/*
 * Copyright (c) 2020-2021 Rob Patro, Avi Srivastava, Hirak Sarkar, Dongze He, Mohsen Zakeri.
 *
 * This file is part of alevin-fry
 * (see https://github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use crate::cellfilter::permit_list_from_file;
use crossbeam_queue::ArrayQueue;
use indicatif::{ProgressBar, ProgressStyle};
#[allow(unused_imports)]
use slog::{crit, info, warn};

use needletail::bitkmer::*;
use sprs::TriMatI;
use std::collections::HashSet;
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use crate::em::{em_optimize_subset, EmInitType};
use crate::utils::read_filter_list;

#[allow(clippy::too_many_arguments)]
pub fn infer(
    //num_bootstraps,
    //init_uniform,
    //summary_stat,
    count_mat_file: String,
    eq_label_file: String,
    usa_mode: bool,
    _use_mtx: bool,
    num_threads: u32,
    filter_list: Option<&str>,
    output_dir: String,
    log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {
    info!(
        log,
        "inferring abundances from equivalence class count input."
    );

    // get the path for the equivalence class count matrix
    let count_mat_path = std::path::Path::new(&count_mat_file);
    let count_mat_parent = count_mat_path
        .parent()
        .unwrap_or_else(|| panic!("cannot get parent path of {:?}", count_mat_path));

    // read the file and convert it to csr (rows are *cells*)
    let count_mat: sprs::CsMatBase<u32, u32, Vec<u32>, Vec<u32>, Vec<u32>, _> =
        match sprs::io::read_matrix_market::<u32, u32, &std::path::Path>(count_mat_path) {
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

    let mut num_cells = count_mat.rows();

    // read in the global equivalence class representation
    let eq_label_path = std::path::Path::new(&eq_label_file);
    let global_eq_classes = Arc::new(crate::eq_class::IndexedEqList::init_from_eqc_file(
        eq_label_path,
    ));

    info!(
        log,
        "read {} equivalence classes from file.",
        global_eq_classes.num_eq_classes()
    );

    // the number of genes (columns) that the output (gene-level) matrix will have
    let num_genes = global_eq_classes.num_genes;

    let usa_offsets = if usa_mode {
        Some(((num_genes / 3) as usize, (2 * num_genes / 3) as usize))
    } else {
        None
    };

    // if we have a filter list, extract it here
    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut retained_bc: HashSet<u64, ahash::RandomState> = HashSet::with_hasher(s);
    let mut filter_bc = false;
    // read the first barcode in the file to get the barcode
    // length.
    let bc_path = count_mat_parent.join("quants_mat_rows.txt");
    let bc_len;
    {
        let mut first_line = String::new();
        let file = match fs::File::open(&bc_path) {
            Ok(file) => file,
            Err(_) => panic!("Unable to read first barcode from {:?}", &bc_path),
        };
        let mut rdr = BufReader::new(file);
        rdr.read_line(&mut first_line).expect("Unable to read line");
        if first_line.ends_with('\n') {
            first_line.pop();
        }
        bc_len = first_line.len() as u16;
    }
    let bc_fname = bc_path
        .to_str()
        .expect("couldn't unwrap barcode file path")
        .to_string();
    let bcvec = permit_list_from_file(bc_fname, bc_len);

    if let Some(fname) = filter_list {
        // read in the fitler list
        match read_filter_list(fname, bc_len) {
            Ok(fset) => {
                // the number of cells we expect to
                // actually process
                num_cells = fset.len();
                retained_bc = fset;
                filter_bc = true;
            }
            Err(e) => {
                return Err(e);
            }
        }
    }

    // the progress bar we'll use to monitor progress of the EM
    let pbar = ProgressBar::new(num_cells as u64);
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
    let cells_to_process = Arc::new(AtomicUsize::new(num_cells));
    // the output cell-by-gene matrix
    let trimat = Arc::new(Mutex::new(TriMatI::<f32, u32>::with_capacity(
        (num_cells, num_genes),
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
                if let Some((cell_num, cell_data)) = in_q.pop() {
                    cells_remaining.fetch_sub(1, Ordering::SeqCst);

                    // given the set of equivalence classes and counts for
                    // this cell (coming from the input matrix), perform
                    // inference to obtain gene-level counts.
                    let counts = em_optimize_subset(
                        &global_eq_classes,
                        &cell_data,
                        &mut unique_evidence,
                        &mut no_ambiguity,
                        EmInitType::Informative,
                        num_genes,
                        false,
                        usa_offsets,
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

    // create our output directory
    let output_path = std::path::Path::new(&output_dir);
    fs::create_dir_all(output_path)?;

    let in_col_path = count_mat_parent.join("quants_mat_cols.txt");
    let out_col_path = output_path.join("quants_mat_cols.txt");
    fs::copy(in_col_path, out_col_path).expect("could not copy column (gene) names to output");

    let out_barcode_path = output_path.join("quants_mat_rows.txt");
    let out_bc_file =
        fs::File::create(out_barcode_path).expect("couldn't create output barcode file");
    let mut bc_writer = BufWriter::new(out_bc_file);

    let mut processed_ind = 0_usize;
    let short_bc_len = bc_len as u8;
    // iterate over the rows (cells) of the equivalence class count
    // matrix.  The iterator will yield a pair of the row index and
    // an iterator over the row data.
    // we zip the bcvec iterator (vector of barcodes in row order)
    // with the actual rows of the matrix.
    for (barcode, row_vec) in bcvec.iter().zip(count_mat.outer_iterator()) {
        let process_cell = if filter_bc {
            // if the reatined_bc list contains this cell id
            // then process it
            retained_bc.contains(barcode)
        } else {
            true
        };
        if process_cell {
            // write to barcode file
            let bc_bytes = &bitmer_to_bytes((*barcode, short_bc_len))[..];
            writeln!(&mut bc_writer, "{}", unsafe {
                std::str::from_utf8_unchecked(bc_bytes)
            })
            .expect("can't write to barcode file.");

            // gather the data from the row's iterator into a vector of
            // (eq_id, count) tuples.
            let cell_data: Vec<(u32, u32)> = row_vec.iter().map(|e| (e.0 as u32, *e.1)).collect();
            // keep pushing this data onto our work queue while we can.
            // launch off these cells on the queue
            let mut cd_clone = (processed_ind, cell_data.clone());
            // keep trying until we can push this payload
            while let Err(t) = q.push(cd_clone) {
                cd_clone = t;
                // no point trying to push if the queue is full
                while q.is_full() {}
            }
            pbar.inc(1);
            processed_ind += 1;
        }
    }

    for h in thread_handles {
        match h.join() {
            Ok(_) => {}
            Err(_e) => {
                info!(log, "thread panicked");
            }
        }
    }

    let pb_msg = format!("finished quantifying {} cells.", num_cells);
    pbar.finish_with_message(pb_msg);

    // finally, we write our output matrix of gene
    // counts.
    let output_matrix_path = output_path.join("quants_mat.mtx");
    let writer_deref = trimat.lock();
    let writer = &*writer_deref.unwrap();
    sprs::io::write_matrix_market(&output_matrix_path, writer)?;

    Ok(())
}
