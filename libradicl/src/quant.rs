// Copyright 2020 Rob Patro, Avi Srivastava. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

extern crate bincode;
extern crate crossbeam_queue;
extern crate fasthash;
extern crate indicatif;
extern crate needletail;
extern crate petgraph;
extern crate serde;
extern crate slog;

//use executors::crossbeam_channel_pool;
//use executors::*;
//use std::sync::mpsc::channel;

use self::indicatif::{ProgressBar, ProgressStyle};
use self::petgraph::prelude::*;
#[allow(unused_imports)]
use self::slog::crit;
use self::slog::info;
use crate as libradicl;
use crossbeam_queue::ArrayQueue;
//use crossbeam_channel::bounded;
use fasthash::sea;
use needletail::bitkmer::*;
use scroll::Pwrite;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io;
use std::io::Read;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use flate2::write::GzEncoder;
use flate2::Compression;

use self::libradicl::em::em_optimize;
use self::libradicl::pugutils;
use self::libradicl::schema::{EqMap, PUGEdgeType, ResolutionStrategy};
use self::libradicl::utils::*;

/// Extracts the parsimonious UMI graphs (PUGs) from the
/// equivalence class map for a given cell.
fn extract_graph(
    eqmap: &EqMap,
    log: &slog::Logger,
) -> petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed> {
    let verbose = false;

    fn get_set() -> HashSet<u32, fasthash::sea::Hash64> {
        HashSet::with_hasher(sea::Hash64)
    }

    // given 2 pairs (UMI, count), determine if an edge exists
    // between them, and if so, what type.
    let has_edge = |x: &(u64, u32), y: &(u64, u32)| -> PUGEdgeType {
        if x.0 == y.0 {
            return PUGEdgeType::BiDirected;
        }
        if x.1 > (2 * y.1 - 1) {
            if count_diff_2_bit_packed(x.0, y.0) < 2 {
                return PUGEdgeType::XToY;
            } else {
                return PUGEdgeType::NoEdge;
            }
        } else if y.1 > (2 * x.1 - 1) {
            if count_diff_2_bit_packed(x.0, y.0) < 2 {
                return PUGEdgeType::YToX;
            } else {
                return PUGEdgeType::NoEdge;
            }
        }
        PUGEdgeType::NoEdge
    };

    let mut _bidirected = 0u64;
    let mut _unidirected = 0u64;

    let mut graph = DiGraphMap::<(u32, u32), ()>::new();

    // for every equivalence class in this cell
    for eqid in 0..eqmap.num_eq_classes() {
        if verbose && eqid % 1000 == 0 {
            print!("\rprocessed {:?} eq classes", eqid);
            io::stdout().flush().expect("Could not flush stdout");
        }
        //ctr += 1;

        // get the info Vec<(UMI, frequency)>
        let eq = &eqmap.eqc_info[eqid];

        // for each (umi, count) pair and its index
        let u1 = &eq.umis;
        for (xi, x) in u1.iter().enumerate() {
            // add a node
            graph.add_node((eqid as u32, xi as u32));

            // for each (umi, freq) pair and node after this one
            for (xi2, x2) in u1.iter().enumerate().skip(xi + 1) {
                //for xi2 in (xi + 1)..u1.len() {
                // x2 is the other (umi, freq) pair
                //let x2 = &u1[xi2];

                // add a node for it
                graph.add_node((eqid as u32, xi2 as u32));

                // determine if an edge exists between x and x2, and if so, what kind
                let et = has_edge(&x, &x2);
                // for each type of edge, add the appropriate edge in the graph
                match et {
                    PUGEdgeType::BiDirected => {
                        graph.add_edge((eqid as u32, xi as u32), (eqid as u32, xi2 as u32), ());
                        graph.add_edge((eqid as u32, xi2 as u32), (eqid as u32, xi as u32), ());
                        _bidirected += 1;
                        //if multi_gene_vec[eqid] == true {
                        //    bidirected_in_multigene += 1;
                        //}
                    }
                    PUGEdgeType::XToY => {
                        graph.add_edge((eqid as u32, xi as u32), (eqid as u32, xi2 as u32), ());
                        _unidirected += 1;
                        //if multi_gene_vec[eqid] == true {
                        //    unidirected_in_multigene += 1;
                        //}
                    }
                    PUGEdgeType::YToX => {
                        graph.add_edge((eqid as u32, xi2 as u32), (eqid as u32, xi as u32), ());
                        _unidirected += 1;
                        //if multi_gene_vec[eqid] == true {
                        //    unidirected_in_multigene += 1;
                        //}
                    }
                    PUGEdgeType::NoEdge => {}
                }
            }
        }

        let mut hset = get_set();

        // for every reference id in this eq class
        for r in eqmap.refs_for_eqc(eqid as u32) {
            // find the equivalence classes sharing this reference
            for eq2id in eqmap.eq_classes_containing(*r).iter() {
                // if eq2id <= eqid, then we already observed the relevant edges
                // when we process eq2id
                if (*eq2id as usize) <= eqid {
                    continue;
                }
                // otherwise, if we have already processed this other equivalence
                // class because it shares _another_ reference (apart from r) with
                // the current equivalence class, then skip it.
                if hset.contains(eq2id) {
                    continue;
                }

                // recall that we processed this eq class as a neighbor of eqid
                hset.insert(*eq2id);
                let eq2 = &eqmap.eqc_info[*eq2id as usize];

                // compare all the umis between eqid and eq2id
                let u2 = &eq2.umis;
                for (xi, x) in u1.iter().enumerate() {
                    // Node for equiv : eqid and umi : xi
                    graph.add_node((eqid as u32, xi as u32));
                    for (yi, y) in u2.iter().enumerate() {
                        // Node for equiv : eq2id and umi : yi
                        graph.add_node((*eq2id as u32, yi as u32));
                        let et = has_edge(&x, &y);
                        match et {
                            PUGEdgeType::BiDirected => {
                                graph.add_edge((eqid as u32, xi as u32), (*eq2id, yi as u32), ());
                                graph.add_edge((*eq2id, yi as u32), (eqid as u32, xi as u32), ());
                                _bidirected += 1;
                                //if multi_gene_vec[eqid] == true
                                //    || multi_gene_vec[*eq2id as usize] == true
                                //{
                                //    bidirected_in_multigene += 1;
                                //}
                            }
                            PUGEdgeType::XToY => {
                                graph.add_edge((eqid as u32, xi as u32), (*eq2id, yi as u32), ());
                                _unidirected += 1;
                                //if multi_gene_vec[eqid] == true
                                //    || multi_gene_vec[*eq2id as usize] == true
                                //{
                                //    unidirected_in_multigene += 1;
                                //}
                            }
                            PUGEdgeType::YToX => {
                                graph.add_edge((*eq2id, yi as u32), (eqid as u32, xi as u32), ());
                                _unidirected += 1;
                                //if multi_gene_vec[eqid] == true
                                //    || multi_gene_vec[*eq2id as usize] == true
                                //{
                                //    unidirected_in_multigene += 1;
                                //}
                            }
                            PUGEdgeType::NoEdge => {}
                        }
                    }
                }
            }
        }
    }

    if verbose {
        info!(
            log,
            "\n\nsize of graph ({:?}, {:?})\n\n",
            graph.node_count(),
            graph.edge_count()
        );
    }
    graph
}

pub fn quantify(
    input_dir: String,
    tg_map: String,
    output_dir: String,
    num_threads: u32,
    resolution: ResolutionStrategy,
    //no_em: bool,
    //naive: bool,
    log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {
    let parent = std::path::Path::new(&input_dir);
    let i_file = File::open(parent.join("map.collated.rad")).unwrap();
    let mut br = BufReader::new(i_file);
    let hdr = libradicl::RADHeader::from_bytes(&mut br);
    info!(
        log,
        "paired : {:?}, ref_count : {:?}, num_chunks : {:?}",
        hdr.is_paired != 0,
        hdr.ref_count,
        hdr.num_chunks
    );

    // now that we have the header, parse and convert the
    // tgmap.

    // first, build a hash of each transcript to it's index
    let mut rname_to_id: HashMap<String, u32> = HashMap::with_capacity(hdr.ref_count as usize);
    for (i, n) in hdr.ref_names.iter().enumerate() {
        rname_to_id.insert(n.clone(), i as u32);
    }
    //println!("{:?}", hdr.ref_names);

    // will hold the unique gene names in the order they are encountered
    let mut gene_names: Vec<String> = Vec::with_capacity((hdr.ref_count / 2) as usize);
    let mut gene_name_to_id: HashMap<String, u32> = HashMap::new();

    // now read in the transcript to gene map
    type TSVRec = (String, String);

    // map each transcript id to the corresponding gene id
    // the transcript name can be looked up from the id in the RAD header,
    // and the gene name can be looked up from the id in the gene_names
    // vector.
    let mut tid_to_gid = vec![u32::MAX; hdr.ref_count as usize];

    let t2g_file = std::fs::File::open(tg_map).expect("couldn't open file");
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(t2g_file);

    // now, map each transcript index to it's corresponding gene index
    let mut found = 0usize;
    for result in rdr.deserialize() {
        let record: TSVRec = result?;
        // first, get the id for this gene
        let next_id = gene_name_to_id.len() as u32;
        let gene_id = *gene_name_to_id.entry(record.1.clone()).or_insert(next_id);
        // if we haven't added this gene name already, then
        // append it now to the list of gene names.
        if gene_id == next_id {
            gene_names.push(record.1);
        }
        // get the transcript id
        if let Some(transcript_id) = rname_to_id.get(&record.0) {
            found += 1;
            tid_to_gid[*transcript_id as usize] = gene_id;
        }
    }
    assert_eq!(
        found, hdr.ref_count as usize,
        "The tg-map must contain a gene mapping for all transcripts in the header"
    );

    info!(
        log,
        "tg-map contained {} genes mapping to {} transcripts.",
        gene_names.len(),
        found
    );

    // file-level
    let fl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    let ft_vals = libradicl::FileTags::from_bytes(&mut br);
    info!(log, "File-level tag values {:?}", ft_vals);

    let bct = rl_tags.tags[0].typeid;
    let umit = rl_tags.tags[1].typeid;

    let mut _num_reads: usize = 0;

    let pbar = ProgressBar::new(hdr.num_chunks);
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );

    // Trying this parallelization strategy to avoid
    // many temporary data structures.

    // We have a work queue that contains collated chunks
    // (i.e. data at the level of each cell).  One thread
    // populates the queue, and the remaining worker threads
    // pop a chunk, perform the quantification, and update the
    // output.  The updating of the output requires acquiring
    // two locks (1) to update the data in the matrix and
    // (2) to write to the barcode file.  We also have to
    // decrement an atomic coutner for the numebr of cells that
    // remain to be processed.

    // create a thread-safe queue based on the number of worker threads
    let n_workers = if num_threads > 1 {
        (num_threads - 1) as usize
    } else {
        1
    };
    let q = Arc::new(ArrayQueue::<(usize, u32, u32, Vec<u8>)>::new(4 * n_workers));

    // the number of cells left to process
    let cells_to_process = Arc::new(AtomicUsize::new(hdr.num_chunks as usize));
    // each thread needs a *read-only* copy of this transcript <-> gene map
    let tid_to_gid_shared = std::sync::Arc::new(tid_to_gid);
    // the number of reference sequences
    let ref_count = hdr.ref_count as u32;
    // the types for the barcodes and umis
    let bc_type = libradicl::decode_int_type_tag(bct).expect("unsupported barcode type id.");
    let umi_type = libradicl::decode_int_type_tag(umit).expect("unsupported umi type id.");
    // the number of genes (different than the number of reference sequences, which are transcripts)
    let num_genes = gene_name_to_id.len();

    // create our output directory
    let output_path = std::path::Path::new(&output_dir);
    fs::create_dir_all(output_path)?;

    // well need a protected handle to write out the barcode
    let bc_path = output_path.join("barcodes.txt");
    let bc_file = fs::File::create(bc_path)?;

    let mat_path = output_path.join("counts.eds.gz");
    let buffered = GzEncoder::new(fs::File::create(mat_path)?, Compression::default());

    let bc_writer = Arc::new(Mutex::new((
        BufWriter::new(bc_file),
        BufWriter::new(buffered),
    )));

    // for each worker, spawn off a thread
    for _worker in 0..n_workers {
        // each thread will need to access the work queue
        let in_q = q.clone();
        // and the logger
        let log = log.clone();
        // the shared tid_to_gid map
        let tid_to_gid = tid_to_gid_shared.clone();
        // and the atomic counter of remaining work
        let cells_remaining = cells_to_process.clone();
        // they will need to know the bc and umi type
        let bc_type = bc_type;
        let umi_type = umi_type;
        // and the file writer
        let bcout = bc_writer.clone();
        // and will need to know the barcode length
        let bclen = ft_vals.bclen;

        // now, make the worker thread
        std::thread::spawn(move || {
            // these can be created once and cleared after processing
            // each cell.
            let mut unique_evidence = vec![false; num_genes];
            let mut no_ambiguity = vec![false; num_genes];
            let mut eq_map = EqMap::new(ref_count);

            // pop from the work queue until everything is
            // processed
            while cells_remaining.load(Ordering::SeqCst) > 0 {
                if let Ok((cell_num, _nbyte, _nrec, buf)) = in_q.pop() {
                    let mut nbr = BufReader::new(&buf[..]);
                    let mut c = libradicl::Chunk::from_bytes(&mut nbr, &bc_type, &umi_type);
                    let bc = c.reads.first().expect("chunk with no reads").bc;
                    eq_map.init_from_chunk(&mut c);

                    let counts: Vec<f32>;
                    match resolution {
                        ResolutionStrategy::CellRangerLike => {
                            counts = pugutils::get_num_molecules_cell_ranger_like(
                                &eq_map,
                                &tid_to_gid,
                                num_genes,
                                &log,
                            );
                        }
                        ResolutionStrategy::Trivial => {
                            counts = pugutils::get_num_molecules_trivial_discard_all_ambig(
                                &eq_map,
                                &tid_to_gid,
                                num_genes,
                                &log,
                            );
                        }
                        ResolutionStrategy::Parsimony => {
                            let g = extract_graph(&eq_map, &log);
                            let gene_eqc =
                                pugutils::get_num_molecules(&g, &eq_map, &tid_to_gid, &log);
                            counts = em_optimize(
                                &gene_eqc,
                                &mut unique_evidence,
                                &mut no_ambiguity,
                                num_genes,
                                true,
                                &log,
                            );
                        }
                        ResolutionStrategy::Full => {
                            let g = extract_graph(&eq_map, &log);
                            let gene_eqc =
                                pugutils::get_num_molecules(&g, &eq_map, &tid_to_gid, &log);
                            counts = em_optimize(
                                &gene_eqc,
                                &mut unique_evidence,
                                &mut no_ambiguity,
                                num_genes,
                                false,
                                &log,
                            );
                        }
                    }
                    // clear our local variables
                    eq_map.clear();
                    // Note: there is a fill method, but it is only on
                    // the nightly branch.  Use this for now:
                    unique_evidence.clear();
                    unique_evidence.resize(num_genes, false);
                    no_ambiguity.clear();
                    no_ambiguity.resize(num_genes, false);
                    // done clearing

                    {
                        // writing the files
                        let bc_mer: BitKmer = (bc, bclen as u8);
                        let eds_bytes: Vec<u8> = sce::eds::as_bytes(counts, num_genes)
                            .expect("can't conver vector to eds");

                        let writer = &mut *bcout.lock().unwrap();
                        // write to barcode file
                        writeln!(&mut writer.0, "{}\t{}", cell_num, unsafe {
                            std::str::from_utf8_unchecked(&bitmer_to_bytes(bc_mer)[..])
                        })
                        .expect("can't write to barcode file.");

                        // write to matrix file
                        writer
                            .1
                            .write_all(&eds_bytes)
                            .expect("can't write to matrix file.");
                    }

                    cells_remaining.fetch_sub(1, Ordering::SeqCst);
                }
            }
        });
    }

    let mut buf = vec![0u8; 65536];
    for cell_num in 0..(hdr.num_chunks as usize) {
        let (nbytes_chunk, nrec_chunk) = libradicl::Chunk::read_header(&mut br);
        buf.resize(nbytes_chunk as usize + 8, 0);
        buf.pwrite::<u32>(nbytes_chunk, 0)?;
        buf.pwrite::<u32>(nrec_chunk, 4)?;
        br.read_exact(&mut buf[8..]).unwrap();
        loop {
            let r = q.push((cell_num, nbytes_chunk, nrec_chunk, buf.clone()));
            if r.is_ok() {
                pbar.inc(1);
                break;
            }
        }
    }

    let gn_path = output_path.join("gene_names.txt");
    let gn_file = File::create(gn_path).expect("couldn't create gene name file.");
    let mut gn_writer = BufWriter::new(gn_file);
    for g in gene_names {
        gn_writer.write_all(format!("{}\n", g).as_bytes())?;
    }

    while cells_to_process.load(Ordering::SeqCst) > 0 {
        // waiting
    }

    Ok(())
}
