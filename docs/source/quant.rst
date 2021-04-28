quant
=====

The ``quant`` command takes a collated RAD file and performs feature (e.g. gene) quantification, outputting
a sparse matrix of de-duplicated counts as well as a list of labels for the rows and columns.  The ``quant``
command takes an input directory containing the collated RAD file, a transcript-to-gene map, an output directory
where the results will be written and a "resolution strategy" (described below).  Quantification is 
multi-threaded, so it also, optionally, takes as an arguments the number of threads to use concurrently.

The transcript-to-gene map should be a two-column (headerless) tab separated file where the first column 
contains a transcript name and the second column contains the corresponding gene name for this transcript.

The ``quant`` command exposes a number of different resolution strategies.  They are:

* ``cr-like`` : This strategy is like the one adopted in cell-ranger, except that it does not first collapse 1-edit-distance UMIs.  Within each cell barcode, a list of (gene, UMI, count) tuples is created. If a read maps to more than one gene, then it generates more than one such tuple.  The tuples are then sorted lexicographically (first by gene id, then by UMI, and then by count).  Any UMI that aligns to only a single gene is assigned to that gene.  UMIs that align to more than one gene are assigned to the gene with the highest count for this UMI.  If there is a tie for the highest count gene for this UMI, then the corresponding reads are simply discarded.

* ``cr-like-em`` : This strategy is like ``cr-like``, except that when a UMI has genes to which it matches with equal frequency, rather than discard the UMIs, the genes are treated as an equivalence class, and the counts for each gene are determined via an expectation maximization algorithm.

* ``full`` : This implements the algorithm described in the alevin_ paper.  Briefly, it builds a graph among the set of reads that align to an overlapping set of transcripts and that have similar (within an edit distance of 1) UMIs.  It then attempts to find a parsimonious cover for this graph using the fewest number of possible transcripts.  If a unique parsimonious cover is found, then the (deduplicated) reads are assigned directly to the genes that yield the most parsimonious cover. If multiple equally-parsimonious covers exist, then the reads are considered multi-mapping at the gene level and they are probabilistically resolved using an expectation maximization (EM) algorithm. 

* ``parsimony`` : This strategy is the same as "full", except that it does *not* probabilistically resolve reads that remain as gene-multimapping after applying the parsimony criterion.  Instead, reads that do not have a unique most-parsimonious assignment are discarded. 

* ``trivial`` : This strategy does not search for 1 edit-distance neighbors of UMIs.  Instead, it first discards any reads that multi-map at the gene level.  The reads that remain then all map uniquely to a single gene.  These reads are deduplicated by (exact) UMI, and the number of distinct UMIs mapping to each gene are taken as that gene's count in the current cell.

Additionally, this command can optionally take the following flags (note that not all resolution strategies are compatible with these flags):

* ``-d, --dump-eqclasses`` : This flag will cause a gene-level, UMI-deduplicated, equivalence class counts file to be written to the output directory in addition to the gene-level count matrix.  This can be used for subsequent analyses where gene-ambiguous reads have been neither resovled nor discarded.

* ``-b, --num-bootstraps`` : This flag will cause bootstrap inferential replicate information to be written to the output directory.  This provides a measure of the inferential uncertainty in the gene-level estimates provided by ``alevin-fry`` when run with a method using the EM algorithm for gene-level abundance estimation.  This information can be used with downstream testing, like differential expression testing using swish.  This flag is only meaningful with the ``cr-like-em`` or ``full`` resolution modes.

* ``--summary-stat`` : This flag will write the summary statistics of the bootstrap replicates (i.e. the mean and variance of the inferential replicates).  This provides the most important information for uncertainty-aware downstream analysis, while requiring much less storage space than the full bootstrap replicate information.  This flag is only meaningful when ``--num-bootstraps`` is meaningful.

* ``--quant-subset <sfile>`` : This optional argument provides a file containing list of barcodes to quantify (one barcode per line, written as a string), those not in this list will be ignored during inference and will not appear in the output quantification matrix.  If this argument is not provided, then all of the original barcodes will be quantified.

* ``--use-mtx`` : This flag will cause the output to be written in matrix market coordinate format rather than in EDS format.

output
------

The output of the ``quant`` command consists of 5 files: ``quants_mat_rows.txt``,
``counts.eds.gz``, ``quants_mat_cols.txt``, ``meta_info.json``, and ``features.txt``. 
The ``meta_info.json`` file contains information about the quantification run,
such as the method used for UMI resolution.  The ``features.txt`` file contains
cell-level information designed to be useful in post-quantification cell filtering
(better determining "true" cells from background, noise, doublets etc.).
The other three files all correspond to quantification information.

The ``counts.eds.gz`` is a file in EDS_ format that stores the gene-by-cell
expression matrix. The two other files provide the labels for the rows and
columns of this matrix. The ``quants_mat_cols.txt`` file is a text file that
contains the names of the rows of the matrix, in the order in which it is
written, with one gene name written per line. The ``quants_mat_rows.txt`` file is a
text file that contains the names of the columns of the matrix, in the order
in which it is written, with one barcode name written per line.

.. _alevin: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y
.. _EDS: https://github.com/COMBINE-lab/EDS

..
  matrix market coordinate format file where the number of *rows* is equal to the number of
  genes and the number of columns is equal to the number of *cells*. The header
  line encodes the number of rows, columns and non-zero entries. The subsequent
  lines (1-based indexing) encode the locations and values of the non-zero
  entries.  This entire ``.mtx`` format file is gzipped during output to minimize
  disk space. 

