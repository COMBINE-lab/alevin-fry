/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use anyhow::{Context, anyhow, bail};
use crossbeam_queue::ArrayQueue;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
#[allow(unused_imports)]
use slog::{crit, info, warn};

use needletail::bitkmer::*;
use sprs::TriMatI;
use std::collections::HashSet;
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use crate::em::{EmInitType, EmSubsetScratch, em_optimize_subset_with_scratch};

/// The cell barcode of a count-matrix row name: `quant` names a row
/// `<sample>_<barcode>` in a multi-sample run, and by its barcode otherwise.
fn row_barcode(name: &str) -> &str {
    name.rsplit_once('_').map_or(name, |(_, bc)| bc)
}

#[allow(clippy::too_many_arguments)]
pub fn infer(
    //num_bootstraps,
    //init_uniform,
    //summary_stat,
    count_mat_file: &PathBuf,
    eq_label_file: &PathBuf,
    usa_mode: bool,
    num_threads: u32,
    filter_list: Option<&PathBuf>,
    output_dir: &PathBuf,
    log: &slog::Logger,
) -> anyhow::Result<()> {
    info!(
        log,
        "inferring abundances from equivalence class count input."
    );

    // get the path for the equivalence class count matrix
    let count_mat_path = std::path::Path::new(count_mat_file);
    let count_mat_parent = count_mat_path
        .parent()
        .unwrap_or_else(|| panic!("cannot get parent path of {:?}", count_mat_path));

    // Read the file and convert it to csr (rows are *cells*).
    //
    // `quant` builds this matrix as `f32`, so the MatrixMarket header it writes
    // says `real`. sprs's reader is strict in *both* directions — a `real` file
    // cannot be read as an integer type and an `integer` file cannot be read as
    // a float one — so asking only for integers made `infer` reject the output
    // of the `quant` it ships alongside (issue #150). Ask for `real` first,
    // then fall back to `integer` so a matrix from an older release or another
    // tool still loads.
    let count_mat: sprs::CsMatI<f32, u32> =
        match sprs::io::read_matrix_market::<f32, u32, &std::path::Path>(count_mat_path) {
            Ok(t) => t.to_csr(),
            Err(float_err) => {
                match sprs::io::read_matrix_market::<i32, u32, &std::path::Path>(count_mat_path) {
                    Ok(t) => t.to_csr().map(|&v| v as f32),
                    Err(int_err) => {
                        warn!(
                            log,
                            "error reading mtx file as real ({:?}) or as integer ({:?})",
                            float_err,
                            int_err
                        );
                        return Err(anyhow!("error reading mtx format matrix : {}", float_err));
                    }
                }
            }
        };

    info!(
        log,
        "read {} x {} equivalence class count matrix.",
        count_mat.rows(),
        count_mat.cols()
    );

    // read in the global equivalence class representation
    let eq_label_path = std::path::Path::new(eq_label_file);
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

    // The output's columns are named by the quants_mat_cols.txt beside the count
    // matrix, so the ids in the equivalence classes must index it. They don't for
    // `quant-probe -d`, whose ids are probes; without this check, `infer` would
    // write a probe-by-cell matrix labelled with gene names.
    let in_col_path = count_mat_parent.join("quants_mat_cols.txt");
    let num_col_names =
        BufReader::new(fs::File::open(&in_col_path).with_context(|| {
            format!("could not open the column names {}", in_col_path.display())
        })?)
        .lines()
        .count();
    if num_col_names != num_genes {
        let probe_col_path = count_mat_parent.join("probe_quants_mat_cols.txt");
        let hint = if probe_col_path.exists() {
            format!(
                " {} is beside it, so these look like the equivalence classes of `quant-probe -d`, whose ids are probes; `infer` produces gene-level counts and does not take them.",
                probe_col_path.display()
            )
        } else {
            String::new()
        };
        bail!(
            "{} has {} features, but {} names {}; the equivalence class ids must index the column names.{}",
            eq_label_path.display(),
            num_genes,
            in_col_path.display(),
            num_col_names,
            hint
        );
    }

    let usa_offsets = if usa_mode {
        Some((num_genes / 3, (2 * num_genes / 3)))
    } else {
        None
    };

    // The row names of the count matrix: a cell barcode, or `<sample>_<barcode>`
    // in a multi-sample run. They are written to the output as they are read.
    let bc_path = count_mat_parent.join("quants_mat_rows.txt");
    let row_names: Vec<String> = BufReader::new(
        fs::File::open(&bc_path)
            .with_context(|| format!("could not open the row names {}", bc_path.display()))?,
    )
    .lines()
    .collect::<Result<_, _>>()
    .with_context(|| format!("could not read the row names {}", bc_path.display()))?;
    if row_names.len() != count_mat.rows() {
        bail!(
            "{} names {} rows, but the count matrix has {}",
            bc_path.display(),
            row_names.len(),
            count_mat.rows()
        );
    }

    // Which rows to quantify. A --quant-subset line that names a row of a
    // multi-sample run (`<sample>_<barcode>`) keeps that row; a bare barcode keeps
    // that barcode in every sample. Bare barcodes are compared as k-mers of the
    // row barcodes' length, as before, so a suffix such as Cell Ranger's `-1`
    // is ignored.
    let retained: Vec<bool> = match filter_list {
        Some(fname) if !row_names.is_empty() => {
            let bc_len = row_barcode(&row_names[0]).len() as u8;
            let kmer = |bc: &str| {
                BitNuclKmer::new(bc.as_bytes(), bc_len, false)
                    .next()
                    .map(|(_, k, _)| k.0)
            };
            let list = fs::read_to_string(fname)
                .with_context(|| format!("could not read {}", fname.display()))?;
            let mut names = HashSet::new();
            let mut barcodes = HashSet::new();
            for line in list.lines().map(str::trim).filter(|l| !l.is_empty()) {
                if line.contains('_') {
                    names.insert(line);
                } else if let Some(k) = kmer(line) {
                    barcodes.insert(k);
                }
            }
            row_names
                .iter()
                .map(|name| {
                    names.contains(name.as_str())
                        || kmer(row_barcode(name)).is_some_and(|k| barcodes.contains(&k))
                })
                .collect()
        }
        _ => vec![true; row_names.len()],
    };
    // The workers stop once this many cells are done, so it must be the number of
    // rows kept, not the length of the --quant-subset list: a listed barcode that
    // no row has would leave them waiting forever.
    let num_cells = retained.iter().filter(|&&keep| keep).count();
    if filter_list.is_some() {
        info!(
            log,
            "quantifying {} of {} cells, from the --quant-subset list",
            num_cells,
            row_names.len()
        );
    }

    // the progress bar we'll use to monitor progress of the EM
    let pbar = ProgressBar::with_draw_target(
        Some(num_cells as u64),
        ProgressDrawTarget::stderr_with_hz(5u8), // update at most 5 times/sec.
    );

    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .expect("ProgressStyle template was invalid.")
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
            let mut em_scratch = EmSubsetScratch::default();

            // Reuse capacity across cells, but grow only to observed
            // expression rather than reserving the full gene axis per worker.
            let mut expressed_vec = Vec::<f32>::new();
            let mut expressed_ind = Vec::<usize>::new();
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
                    let counts = em_optimize_subset_with_scratch(
                        &global_eq_classes,
                        &cell_data,
                        &mut unique_evidence,
                        &mut no_ambiguity,
                        EmInitType::Informative,
                        num_genes,
                        false,
                        usa_offsets,
                        &log,
                        &mut em_scratch,
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
                            writer.add_triplet(row_index, *ind, *val);
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
    let output_path = std::path::Path::new(output_dir);
    fs::create_dir_all(output_path)?;

    let out_col_path = output_path.join("quants_mat_cols.txt");
    fs::copy(in_col_path, out_col_path).expect("could not copy column (gene) names to output");

    let out_barcode_path = output_path.join("quants_mat_rows.txt");
    let out_bc_file =
        fs::File::create(out_barcode_path).expect("couldn't create output barcode file");
    let mut bc_writer = BufWriter::new(out_bc_file);

    let mut processed_ind = 0_usize;
    // iterate over the rows (cells) of the equivalence class count
    // matrix, with their names and whether each is kept.
    for ((name, keep), row_vec) in row_names
        .iter()
        .zip(retained.iter())
        .zip(count_mat.outer_iterator())
    {
        if *keep {
            writeln!(&mut bc_writer, "{}", name).expect("can't write to barcode file.");

            // gather the data from the row's iterator into a vector of
            // (eq_id, count) tuples.
            // NOTE: the matrix market format (and sprs' implementation of it)
            // does not carry an unsigned integer type, so the matrix is read as
            // `f32` and the counts are converted here. These are whole UMI
            // counts held in a float, so round rather than truncate: a value
            // that round-tripped as 4.999999 should become 5, not 4.
            let cell_data: Vec<(u32, u32)> = row_vec
                .iter()
                .map(|e| (e.0 as u32, e.1.round() as u32))
                .collect();
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
    sprs::io::write_matrix_market(output_matrix_path, writer)?;

    Ok(())
}
