collate
=======

This command takes as input a directory containing a RAD file (created by running alevin with the ``--justAlign`` and/or ``--sketch`` flags), as well as the directory generated as the result of running the ``generate-permit-list`` command of ``alevin-fry``, and it will produce an output RAD file that is *collated* by (corrected) cellular barcode.  The collated RAD file can then be quantified with the ``alevin-fry`` ``quant`` command.  It also takes two other arguments (described below) that dictate how the collation and filtering will be performed.

* ``-r, --rad-dir <rad-dir>`` : The directory containing the RAD file to be collated.  This is the *same* directory on which you have previously run ``generate-permit-list`` and that was obtained by running ``alevin`` with the ``--justAlign`` flag).

* ``-i, --input-dir <input-dir>`` : The input directory.  This is the directory that was the *output* of ``generate-permit-list``.  This directory contains information computed by the ``generate-permit-list`` command that will allow successful collation and barcode correction.  This is also the directory where the collated RAD file will be *output*.

* ``--compress`` : This optional flag will tell ``alevin-fry`` to compress the output collated RAD file.  The file will be compressed using the `Snappy compression format <https://github.com/google/snappy/blob/master/format_description.txt>`__ (via the excellent `snap <https://docs.rs/snap/>`__ crate.  If this option is passed, the output file will be written to ``map.collated.rad.sz`` rather than ``map.collated.rad``, and the corresponding status of the file's compression will be written to ``collate.json`` in the output file.  *Note*: The choice to use compression or not has no effect on the final result or the correctness of the output, but it may have some moderate performance implications.  Specifically, it is potentially worth using this flag if you want to minimize disk space, and if you are using a sufficiently large number of threads (as compression happens in parallel, a sufficient number of threads will allow the compressed RAD file to be generated as quickly as the uncompressed).  However, because some internal buffers must be duplicated during parallel compression, the collate step can use a bit more memory if run with the ``--compress`` flag, though the memory usage should still be small and stable over different sized inputs.  There can also be an effect on quantification speed (since the collated RAD file will be decompressed on the fly during quantification), but it should be small since Snappy decompresses very fast, and decompression will only be the limiting factor if you are using a simple resolution strategy (e.g. naive or cr-like) and many quantification threads.
 
* ``--memory-limit <size>`` : The collation buffer budget, excluding correction indexes and the operating-system page cache.  Values such as ``2GiB`` and ``4GB`` are accepted; the default is 2 GiB.  ``--max-records`` remains as a deprecated approximate compatibility control.

Current collation makes one scatter pass over the input RAD file.  GPL's
compiled correction decisions are fused with output bucket identifiers, so a
record is corrected and routed with one direct lookup.  Worker-local temporary
spools are then gathered into the output.  The number of physical spools is
bounded by the worker count rather than by the number of logical cell buckets,
and the correction index is released before the memory-intensive gather phase.

For multi-barcode data, the resulting RAD is hierarchically ordered: records
for each sample are contiguous, and all records for a corrected cell are
collated within that sample.  The manifest uses a complete ordering including
sample and cell-barcode tie-breaks.

If ``correction_plan.bin`` is absent, collation uses the historical correction
files and emits a warning.  A present but malformed or unsupported plan is an
error; it is never silently ignored.

output
------

The ``collate`` command will output all files it creates in the expected format in the output directory that is specified. It will write a file name ``map.collated.rad`` (or ``map.collated.rad.sz`` if run with the ``--compress`` flag), one named ``unmapped_bc_count_collated.bin``, and one named ``collate.json`` in the directory specified by ``-i``.
