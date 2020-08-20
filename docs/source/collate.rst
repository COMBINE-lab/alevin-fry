collate
=======

This command takes as input a RAD file (created by running alevin with the ``--justAlign`` flag), as well as a 
directory generated as the result of running the ``generate-permit-list`` command of ``alevin-fry``, and it will
produce an output RAD file that is *collated* by (corrected) cellular barcode.  The collated RAD file can then 
be quantified with the ``alevin-fry`` ``quant`` command.  It also takes two other arguments (described below) that 
dictate how the collation and filtering will be performed.

* ``-r, --rad-file <rad-file>`` : The RAD file to be collated.  This is the *same* file on which you have previously run ``generate-permit-list`` and that was obtained by running ``alevin`` with the ``--justAlign`` flag).

* ``-i, --input-dir <input-dir>`` : The input directory.  This is the directory that was the *output* of ``generate-permit-list``.  This directory contains information computed by the ``generate-permit-list`` command that will allow successful collation and barcode correction.  This is also the directory where the collated RAD file will be *output*.

* ``-m, --max-records <max-records>`` : The maximum number of read records to keep in memory at once during collation. The ``collate`` command will pass over the input RAD file multiple times collecting the records associated with a set of (corrected) cellular barcodes so that they can be written out in collated format to the output RAD file.  This parameter determines (approximately) how many records will be held in memory at once, and therefore determines the memory usage of the ``collate`` command.  The larger the value used the faster the collation process will be, since fewer passes are made.  The smaller this value, the lower the memory usage will be, at the cost of more passes.  The default value is 10,000,000.  Note that this determines the number of records *approximately*, because a specific barcode will never be split across multiple collation passes.  The algorithm employed is to collect the reads associated with different cellular barcodes in the current pass until the number of reads to be collected *first exceeds* this value.

output
------

The ``collate`` command only has a single output.  It will write a file name
``map.collated.rad`` in the directory specified by ``-i``.