infer
=====

The ``infer`` command takes a gene-resolved equivalence class count matrix, along with the 
equivalence class description file (both output when running ``quant`` with the ``-d`` flag),
and performs inference via the EM to reduce the matrix to a cell-by-gene level count matrix.

This functionality makes it possible to separate the UMI resolution step (which may result in 
UMIs being assigned to equivalence classes of genes rather than individual genes), and the gene-level 
estimation step, that attempts to resolve gene-multimapping UMIs.

T command can takes the following options :

* ``-c, --count-mat <eqc-mat>`` : This provides the path to the (mtx format) matrix of cells by equivalence class counts.

* ``-e, --eq-labels <eq-labels>`` : This provides the path to the file containing the gene labels of the equivalence class description.

* ``-o, --output-dir <output-dir>`` : This provides the output file where the quantification matrix will be written.

* ``-t, --threads <threads>`` : This option provides the number of threads to use for processing [default: number of hardware threads].

output
------

The output of the ``infer`` command consists of a single file; the cell-by-gene count matrix derived from the input 
cell-by-equivalence-class count matrix.  Crucially, the ordering of the rows in the cell-by-gene matrix is the same 
as that in the input matrix, so that another "features" file is not written.
