extern crate slog;
extern crate serde;
extern crate sce;
extern crate sprs;
extern crate fasthash;

#[allow(unused_imports)]
use fasthash::sea::Hash64;
use serde::{Deserialize, Serialize};
use self::slog::{crit, info};
use std::fs::File;
use std::io::{prelude::*, BufReader, BufWriter, Write};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::str;
use std::collections::HashMap;

use crate as libradicl;
use self::libradicl::em::{em_optimize, run_bootstrap};


#[derive(Serialize, Deserialize)]
struct Group {
    ids: Vec<String>,
    groups: Vec<Vec<String>>,
}

#[allow(dead_code)]
type BufferedGZFile = BufWriter<GzEncoder<File>>;
#[allow(dead_code)]
struct BootstrapHelper {
    bsfile: Option<BufferedGZFile>,
    mean_var_files: Option<(BufferedGZFile, BufferedGZFile)>,
}
#[allow(dead_code)]
impl BootstrapHelper {
    fn new(
        output_path: &std::path::Path,
        num_bootstraps: u32,
        summary_stat: bool,
    ) -> BootstrapHelper {
        if num_bootstraps > 0 {
            if summary_stat {
                let bootstrap_mean_path = output_path.join("bootstraps_mean.eds.gz");
                let bootstrap_var_path = output_path.join("bootstraps_var.eds.gz");
                let bt_mean_buffered = GzEncoder::new(
                    File::create(bootstrap_mean_path).unwrap(),
                    Compression::default(),
                );
                let bt_var_buffered = GzEncoder::new(
                    File::create(bootstrap_var_path).unwrap(),
                    Compression::default(),
                );
                BootstrapHelper {
                    bsfile: None,
                    mean_var_files: Some((
                        BufWriter::new(bt_mean_buffered),
                        BufWriter::new(bt_var_buffered),
                    )),
                }
            } else {
                let bootstrap_path = output_path.join("bootstraps.eds.gz");
                let bt_buffered = GzEncoder::new(
                    File::create(bootstrap_path).unwrap(),
                    Compression::default(),
                );
                BootstrapHelper {
                    bsfile: Some(BufWriter::new(bt_buffered)),
                    mean_var_files: None,
                }
            }
        } else {
            BootstrapHelper {
                bsfile: None,
                mean_var_files: None,
            }
        }
    }
}


pub fn calc_diff(
    input_dir: String,
    _num_threads: u32,
    _num_bootstraps : u32,
    _summary_stat: bool,
    group_file_name: String,
    log: &slog::Logger,
){

    // reading the group file
    let reader = BufReader::new(File::open(group_file_name)
        .expect("error")); 
    let groups : Group = serde_json::from_reader(reader).unwrap();
  
    assert_eq!(groups.ids.len(),groups.groups.len());

    // reading matrix file
    let input_path = std::path::Path::new(&input_dir);
    let input_matrix_path = input_path.join("alevin");

    let bc_path = input_matrix_path.join("quants_mat_rows.txt");
    let gn_path = input_matrix_path.join("quants_mat_cols.txt");
    let mat_path = input_matrix_path.join("quants_mat.gz");

    // let boot_helper = BootstrapHelper::new(input_path, num_bootstraps, summary_stat);

    let mut bc_list = Vec::<String>::new();
    let mut gn_list = Vec::<String>::new();
    let mut cell_id_map : HashMap<String,usize> = HashMap::new();
    {
        let bc_buf = BufReader::new(File::open(bc_path).expect("error"));
        for (i,line) in bc_buf.lines()
                       .map(|l| l.expect("Could not parse line"))
                       .enumerate()
        {
            bc_list.push(line.clone());
            cell_id_map.insert(line.clone(), i);
        }
        
        let gn_buf = BufReader::new(File::open(gn_path).expect("error"));
        for line in gn_buf.lines()
                       .map(|l| l.expect("Could not parse line"))
        {
            gn_list.push(line);
        }
    }
    let sce_eds = sce::SingleCellExperiment::from_eds(
        &mat_path.to_str().unwrap(),
        bc_list.clone(),
        gn_list.clone(),
    ).unwrap();

    info!(log,
        "Quant matrix {:?}",
        sce_eds.shape(),
    );

    // read gene equivalence classes
    let mtx_path = input_path.join("geqc_counts.mtx");
    let eq_mat: sprs::TriMatI<f32, u32> =  sprs::io::read_matrix_market(&mtx_path)
                                            .expect("could not write geqc_counts.mtx");
    
    assert_eq!(eq_mat.rows(), bc_list.len());
    let num_eq_classes = eq_mat.cols();
    let num_genes = gn_list.clone().len();
    // let num_cells = bc_list.clone().len();

    info!(
        log,
        "Number of equivalence classes {}",
        num_eq_classes,
    );
    let mut global_eqc : Vec<Vec<u32>> = vec![Vec::new(); num_eq_classes];
    let eqc_path = input_path.join("gene_eqclass.txt.gz");
    {
        let gn_eq_buf = BufReader::new(
            GzDecoder::new(File::open(eqc_path).expect("error"))
        );
        let mut first_line : bool = true;
        for line in gn_eq_buf.lines()
            .map(|l| l.expect("Could not parse line"))
        {
            if first_line {
                first_line = false;
                continue;
            }
            let mut tokens : Vec<u64> = line.split("\t")
                .map(|x| x.parse::<u64>().unwrap()).collect();
            let eq_id : usize = tokens.pop().unwrap() as usize;
            let tokens_u32 : Vec<u32> = tokens.iter().map(|&e| e as u32).collect();
            global_eqc[eq_id] = tokens_u32.clone();
        }
    }

    // Segementing the global eqclass vec
    // make a vector of global vectors
    //let mut global_eq_vec = vec![Vec::new(); groups.ids.len()];
    let num_groups = groups.ids.len();
    let mut global_eq_hash : HashMap<usize,Vec<f32>> = HashMap::new();
    {
        let eq_csc = eq_mat.to_csc();
        let mut cell_to_group : HashMap<usize,usize> = HashMap::new();
        for (group_id, cell_vec) in groups.groups.iter().enumerate() {
            for cell_name in cell_vec.iter(){
                let cell_id = cell_id_map.get(cell_name).unwrap();
                cell_to_group.insert(*cell_id, group_id);
            }
        }

        // let mut eq_vec = vec![0f32;num_eq_classes];
        for (row_ind, row) in eq_csc.outer_iterator().enumerate(){
            // Which group this cell belongs belongs to?
            if let Some(group_id) = cell_to_group.get(&row_ind) {
                if let Some(mut eq_vec) = global_eq_hash.get_mut(&group_id){
                    for (col_ind, x) in row.iter() {
                        // merge the count of all cells
                        eq_vec[col_ind] += *x ;
                    }
                }else{
                    let mut eq_vec = vec![0f32; num_eq_classes];
                    for (col_ind, x) in row.iter() {
                        // merge the count of all cells
                        eq_vec[col_ind] = *x ;
                    }
                    global_eq_hash.insert(*group_id, eq_vec.clone());
                }
            }
        }
    }

    {
        // For each group we need to re-run
        let mat_path = input_path.join("group_mat.gz");
        let mut group_eds_file = BufWriter::new(
            GzEncoder::new(
                File::create(&mat_path).unwrap(), 
                Compression::default()
            )
        );
        for group_id in 0..num_groups {
            // let group_name = groups.ids[group_id].clone();

            // let boot_helper = BootstrapHelper::new(output_path, num_bootstraps, summary_stat);
            // let buffered = GzEncoder::new(fs::File::create(&mat_path)?, Compression::default());

            let eq_vec = global_eq_hash.get(&group_id).unwrap();
            // re-run the EM algorithm for a group of cells on the joint equivalence classes
            // create the hashmap so that we can call em::em_optimize
            let s = fasthash::RandomState::<Hash64>::new();
            let mut gene_eqc: HashMap<Vec<u32>, u32, fasthash::RandomState<Hash64>> =
                HashMap::with_hasher(s);
            for (eq_id, umi_count) in eq_vec.iter().enumerate(){
                gene_eqc.insert(global_eqc[eq_id].clone(), (*umi_count) as u32);
            }
            // run em on the merged equivalence classes
            let counts : Vec<f32>;
            let mut unique_evidence = vec![false; num_genes];
            let mut no_ambiguity = vec![false; num_genes];
            counts = em_optimize(
                &gene_eqc,
                &mut unique_evidence,
                &mut no_ambiguity,
                num_genes,
                true,
                &log,
            );
            let eds_bytes = sce::eds::as_bytes(&counts, num_genes)
                .expect("can't convert vector to eds");
            group_eds_file.write_all(&eds_bytes).expect("can't write matrix file");
        }
    }
}