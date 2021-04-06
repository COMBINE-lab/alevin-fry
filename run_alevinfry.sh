#!/bin/bash

# build expanded txome
# expanded.tx2gene.tsv, expanded.fa and expanded.gtf will be generated and used later
Rscript spliced_unspliced_txome_builder.R \
--fa=path_to_genome_file \
--gtf=path_to_gtf_file \
--outdir= path_to_output_directory

# build salmon index
# expanded.fa is used
path_to_salmon-1.4.0 index \
-t path_to_expanded_fa \
-i path_to_expected_output_directory \
--gencode -p num_threads

# generate map.rad file
# expanded.tx2gene.tsv is used
path_to_salmon-1.4.0 alevin -l ISR \
-i path_to_salmon_index_folder \
-1 path_to_fastq_1 \
-2 path_to_fastq_2 \
-o path_to_expected_outdir \
--tgMap path_to_expanded_tx2gene_file \
-p num_threads --chromium \
--justAlign --sketchMode

# cell filtering
# -v is the velocity flag
path_to_alevinfry generate-permit-list \
-i path_to_rad_file \
-o path_to_expected_fry_outdir \
-k -d fw -v

# collate cells
path_to_alevinfry collate \
-r path_to_rad_file \
-i path_to_fry_outdir \
-t numthreads

# quantification
# -d is the dump_eqclasses flag
# expanded.tx2gene.tsv is used
path_to_alevinfry quant \
-i path_to_collate_outdir \
-o path_to_fry_outdir \
-m path_to_expanded_tx2gene_file \
-t num_threads --use-mtx -d -r parsimony
