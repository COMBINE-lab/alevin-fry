// Copyright 2020 Rob Patro, Hirak Sarkar, Avi Srivastava. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

extern crate ahash;
extern crate bincode;
extern crate crossbeam_queue;
extern crate fasthash;
extern crate indicatif;
extern crate serde;
extern crate slog;

use self::indicatif::{ProgressBar, ProgressStyle};
#[allow(unused_imports)]
use self::slog::{crit, info, warn};
use crate as libradicl;
use crossbeam_queue::ArrayQueue;

// use fasthash::sea;
use sprs::TriMatI;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use self::libradicl::em::em_optimize_subset;

pub fn infer(
    //num_bootstraps,
    //init_uniform,
    //summary_stat,
    count_mat_file: String,
    eq_label_file: String,
    _use_mtx: bool,
    num_threads: u32,
    output_dir: String,
    log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {
    info!(
        log,
        "inferring abundances from equivalence class count input."
    );

    // get the path for the equivalence class count matrix
    let count_mat_path = std::path::Path::new(&count_mat_file);
    // read the file and convert it to csr (rows are *cells*)
    let count_mat = match sprs::io::read_matrix_market::<u32, u32, &std::path::Path>(count_mat_path)
    {
        Ok(t) => t.to_csr(),
        Err(e) => {
            warn!(log, "error reading mtx file{:?}", e);
            return Err(Box::new(e));
        }
    };

    info!(
        log,
        "read {} x {} equivalence class count matrix.",
        count_mat.rows(),
        count_mat.cols()
    );

    // read in the global equivalence class representation
    let eq_label_path = std::path::Path::new(&eq_label_file);
    let global_eq_classes = Arc::new(libradicl::schema::IndexedEqList::init_from_eqc_file(
        eq_label_path,
    ));

    info!(
        log,
        "read {} equivalence classes from file.",
        global_eq_classes.num_eq_classes()
    );

    // the number of genes (columns) that the output (gene-level) matrix will have
    let num_genes = global_eq_classes.num_genes;

    // the progress bar we'll use to monitor progress of the EM
    let pbar = ProgressBar::new(count_mat.rows() as u64);
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );

    // create a thread-safe queue based on the number of worker threads
    let n_workers = if num_threads > 1 {
        (num_threads - 1) as usize
    } else {
        1
    };
    // the queue will hold tuples of the
    // cell id (so that the output matrix is in the same order as input)
    // vector of eq_id and count for each cell
    let q = Arc::new(ArrayQueue::<(usize, Vec<(u32, u32)>)>::new(4 * n_workers));

    let mut thread_handles: Vec<thread::JoinHandle<_>> = Vec::with_capacity(n_workers);

    // the number of cells left to process
    let cells_to_process = Arc::new(AtomicUsize::new(count_mat.rows()));
    // the output cell-by-gene matrix
    let trimat = Arc::new(Mutex::new(TriMatI::<f32, u32>::with_capacity(
        (count_mat.rows(), num_genes),
        count_mat.nnz(),
    )));

    // for each worker, spawn off a thread
    for _worker in 0..n_workers {
        // each thread will need to access the work queue
        let in_q = q.clone();
        // and the logger
        let log = log.clone();
        // and the atomic counter of remaining work
        let cells_remaining = cells_to_process.clone();
        // and the output matrix
        let matout = trimat.clone();
        // and the global set of eq class labels
        let global_eq_classes = global_eq_classes.clone();

        //let unmapped_count = bc_unmapped_map.clone();
        //let mmrate = mmrate.clone();

        // now, make the worker thread
        let handle = std::thread::spawn(move || {
            // these can be created once and cleared after processing
            // each cell.
            let mut unique_evidence = vec![false; num_genes];
            let mut no_ambiguity = vec![false; num_genes];

            let mut expressed_vec = Vec::<f32>::with_capacity(num_genes);
            let mut expressed_ind = Vec::<usize>::with_capacity(num_genes);
            //let mut eds_bytes = Vec::<u8>::new();
            //let mut bt_eds_bytes: Vec<u8> = Vec::new();
            //let mut eds_mean_bytes: Vec<u8> = Vec::new();
            //let mut eds_var_bytes: Vec<u8> = Vec::new();

            // pop from the work queue until everything is
            // processed
            while cells_remaining.load(Ordering::SeqCst) > 0 {
                if let Ok((cell_num, cell_data)) = in_q.pop() {
                    cells_remaining.fetch_sub(1, Ordering::SeqCst);

                    // given the set of equivalence classes and counts for
                    // this cell (coming from the input matrix), perform
                    // inference to obtain gene-level counts.
                    let counts = em_optimize_subset(
                        &global_eq_classes,
                        &cell_data,
                        &mut unique_evidence,
                        &mut no_ambiguity,
                        num_genes,
                        false,
                        &log,
                    );

                    // Note: there is a fill method, but it is only on
                    // the nightly branch.  Use this for now:
                    unique_evidence.clear();
                    unique_evidence.resize(num_genes, false);
                    no_ambiguity.clear();
                    no_ambiguity.resize(num_genes, false);
                    // done clearing

                    let mut max_umi = 0.0f32;
                    let mut _sum_umi = 0.0f32;
                    let mut _num_expr: u32 = 0;
                    expressed_vec.clear();
                    expressed_ind.clear();
                    // go over the output counts and fill in the vector of
                    // expressed gene ids and their corresponding counts
                    for (gn, c) in counts.iter().enumerate() {
                        max_umi = if *c > max_umi { *c } else { max_umi };
                        _sum_umi += *c;
                        if *c > 0.0 {
                            _num_expr += 1;
                            expressed_vec.push(*c);
                            expressed_ind.push(gn);
                        }
                    }

                    /*
                    let num_mapped = nrec;
                    let dedup_rate = sum_umi / num_mapped as f32;

                    let num_unmapped = match unmapped_count.get(&bc) {
                        Some(nu) => *nu,
                        None => 0u32,
                    };

                    let mapping_rate = num_mapped as f32 / (num_mapped + num_unmapped) as f32;

                    // mean of the "expressed" genes
                    let mean_expr = sum_umi / num_expr as f32;
                    // number of genes with expression > expressed mean
                    let num_genes_over_mean = expressed_vec.iter().fold(0u32, |acc, x| {
                        if x > &mean_expr {
                            acc + 1u32
                        } else {
                            acc
                        }
                    });
                    // expressed mean / max expression
                    let mean_by_max = mean_expr / max_umi;
                    */

                    let row_index: usize; // the index for this row (cell)
                    {
                        row_index = cell_num;
                        let writer_deref = matout.lock();
                        let writer = &mut *writer_deref.unwrap();

                        // fill out the triplet matrix in memory
                        for (ind, val) in expressed_ind.iter().zip(expressed_vec.iter()) {
                            writer.add_triplet(row_index as usize, *ind, *val);
                        }
                        /*
                        writeln!(
                            &mut writer.feature_file,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            unsafe { std::str::from_utf8_unchecked(&bc_bytes) },
                            (num_mapped + num_unmapped),
                            num_mapped,
                            sum_umi,
                            mapping_rate,
                            dedup_rate,
                            mean_by_max,
                            num_expr,
                            num_genes_over_mean
                        )
                        .expect("can't write to feature file");
                            */
                        /*
                        if num_bootstraps > 0 {
                            if summary_stat {
                                if let Some((meanf, varf)) =
                                    &mut writer.bootstrap_helper.mean_var_files
                                {
                                    meanf
                                        .write_all(&eds_mean_bytes)
                                        .expect("can't write to bootstrap mean file.");
                                    varf.write_all(&eds_var_bytes)
                                        .expect("can't write to bootstrap var file.");
                                }
                            } else if let Some(bsfile) = &mut writer.bootstrap_helper.bsfile {
                                bsfile
                                    .write_all(&bt_eds_bytes)
                                    .expect("can't write to bootstrap file");
                            }
                        } // done bootstrap writing
                        */
                    }
                } // while we can get work
            } // while cells remain
        });
        thread_handles.push(handle);
    }

    // iterate over the rows (cells) of the equivalence class count
    // matrix.  The iterator will yield a pair of the row index and
    // an iterator over the row data.
    for (row_ind, row_vec) in count_mat.outer_iterator().enumerate() {
        // gather the data from the row's iterator into a vector of
        // (eq_id, count) tuples.
        let cell_data: Vec<(u32, u32)> = row_vec.iter().map(|e| (e.0 as u32, *e.1)).collect();
        // keep pushing this data onto our work queue while we can.
        loop {
            if !q.is_full() {
                let r = q.push((row_ind, cell_data.clone()));
                // if we were successful in pushing, then increment
                // the progress bar.
                if r.is_ok() {
                    pbar.inc(1);
                    break;
                }
            }
        }
    }

    // finally, we write our output matrix of gene
    // counts.
    let writer_deref = trimat.lock();
    let writer = &*writer_deref.unwrap();
    let output_path = std::path::Path::new(&output_dir);
    sprs::io::write_matrix_market(&output_path, writer)?;

    Ok(())
}
