quant
=====

The ``quant'' command takes a collated RAD file and performs feature (e.g. gene) quantification, outputting
a sparse matrix of de-duplicated counts as well as a list of labels for the rows and columns.  The ``quant``
command takes an input directory containing the collated RAD file, a transcript-to-gene map, an output directory
where the results will be written and a "resolution strategy" (described below).  Quantification is 
multi-threaded, so it also, optionally, takes as an arguments the number of threads to use concurrently.

The transcript-to-gene map should be a two-column (headerless) tab separated file where the first column 
contains a transcript name and the second column contains the corresponding gene name for this transcript.

The ``quant`` command exposes a number of different resolution strategies.  They are:

* ``full`` : This is the default resolution strategy.  It implements the algorithm described in the alevin_ 
paper.  Briefly, it builds a graph among the set of reads that align to an overlapping set of transcripts 
and that have similar (within an edit distance of 1) UMIs.  It then attempts to find a parsimonious cover 
for this graph using the fewest number of possible transcripts.  If a unique parsimonious cover is found,
then the (deduplicated) reads are assigned directly to the genes that yield the most parsimonious cover.
If multiple equally-parsimonious covers exist, then the reads are considered multi-mapping at the gene 
level and they are probabilistically resolved using an expectation maximization (EM) algorithm.

* ``parsimony`` : This strategy is the same as ``full'', except that it does *not* probabilistically resolve
reads that remain as gene-multimapping after applying the parsimony criterion.  Instead, reads that do 
not have a unique most-parsimonious assignment are discarded.

* ``trivial`` : This strategy does not search for 1 edit-distance neighbors of UMIs.  Instead, it first 
discards any reads that multi-map at the gene level.  The reads that remain then all map uniquely to a 
single gene.  These reads are deduplicated by (exact) UMI, and the number of distinct UMIs mapping to 
each gene are taken as that gene's count in the current cell.

* ``cr-like`` : This strategy is like the one adopted in cell-ranger, except that it does not first 
collapse 1-edit-distance UMIs.  Within each cell barcode, a list of (gene, UMI, count) tuples is created.
If a read maps to more than one gene, then it generates more than one such tuple.  The tuples are then 
sorted lexicographically (first by gene id, then by UMI, and then by count).  Any UMI that aligns to only 
a single gene is assigned to that gene.  UMIs that align to more than one gene are assigned to the gene 
with the highest count for this UMI.  If there is a tie for the highest count gene for this UMI, then the 
corresponding reads are simply discarded.

output
------

The output of the ``quant`` command consists of 3 files: ``barcodes.txt``,
``counts.mtx`` and ``gene_names.txt``. The ``counts.mtx`` is a matrix market
coordinate format file where the number of *rows* is equal to the number of
genes and the number of columns is equal to the number of *cells*. The header
line encodes the number of rows, columns and non-zero entries. The subsequent
lines (1-based indexing) encode the locations and values of the non-zero
entries. The two other files provide the labels for the rows and columns of
this matrix. The ``gene_names.txt`` file is a text file that contains the
names of the rows of the matrix, in the order in which it is written, with
one gene name written per line. The ``barcodes.txt`` file is a text file that
contains the names of the columns of the matrix, in the order in which it is
written, with one barcode name written per line.



.. _alevin: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y