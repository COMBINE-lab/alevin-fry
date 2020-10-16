generate-permit-list
====================

This command takes as input an output directory containing a RAD file (created by running alevin with the ``--justAlign`` flag), 
and it determines what cell barcodes should be associated with "true" cells, which should be corrected to
some "true" barcode, and which should simply be ignored / discarded. This
command has 4 required arguments; the path to an input directory ``--input``,
the path to an output directory ``--output-dir`` (which will be created if it
doesn't exist), the expected orientation of properly mapped reads
``--expected-ori`` (the options are 'fw' (filters out alignments to the
reverse complement strand), 'rc' (filter out alignments to the forward
strand) and 'both' or 'either' (do not filter any alignments)), and then one
of the following mutually exclusive options (which determines how the "true"
barcodes are decided):

* ``--knee-distance``: This flag will use the distance method that is used in the whitelist command of 
  UMI-tools to attempt to automatically determine the number of true barcodes. Briefly, this 
  method first counts the number of reads associated with each barcode, and then sorts the barcodes in 
  descending order by their associated read count. It then constructs the cumulative distribution function 
  from this sorted list of frequencies. Finally, it applies an iterative algorithm to attempt to determine the optimal 
  number of barcodes to include by looking for a "knee" or "elbow" in the CDF graph. The algorithm considers 
  each barcode in the CDF where it's x-coordinate is equal to this barcode's rank divided by the total number 
  of barcodes (i.e. its normalized rank) and the y-coordinate is equal to the (normalized) cumulative frequency achieved 
  at this barcode. It then computes the distance of this barcode from the line x=y 
  (defined by the start and end of the CDF). The initial knee is predicted as the point that has the maximum distance 
  from the x=y line. The algorithm is iterative, because experiments with many low-quality barcodes may predict too many 
  valid barcodes using this method. Thus, the algorithm is run repeatedly, each time considering a prefix of the CDF from 
  index 0 through the previous knee's index * 5. Once two subsequent iterations of the algorithm return the same 
  knee point, the algorithm terminates.

* ``--force-cells <ncells>``: This option will count the number of reads associated with each barcode, and sort the barcodes 
  in descending order of frequency. Then, it will consider the first <ncells> barcodes to be valid. Any barcode that has 
  a number of reads >= to the <ncells>-th barcode will be considered part of the permit list, all others will not 
  (but will be considered for correction to this permit list).

* ``--valid-bc <bcfile>``: This option will read the provided file <bcfile> and treat it as an explicitly-provided list of true 
  barcodes. Barcodes appearing in this list will be considered true, and barcodes will be corrected to this list.

* ``--expect-cells <ncells>``: Not currently implemented.

output
------

The ``generate-permit-list`` command outputs a number of different files in the output directory.  Not all files are 
relevant to users of ``alevin-fry``, but the files are described here.

1. The file ``all_freq.tsv`` is a two-column tab-separated file that lists, for each distinct barcode in the input RAD file, the number of read records that were tagged with this barcode.

2. The file ``permit_freq.tsv`` is a two-column tab-separated file that lists, for each barcode in the input RAD file that is determined to be a *true* barcode, the number of read records associated with this barcode.

3. The file ``permit_map.bin`` is a binary file (a serde serialized HashMap) that maps each barcode in the input RAD file that is within an edit distance of 1 to some *true* barcode to the barcode to which it corrects.  This allows the ``collate`` command to group together all of the read records corresponding to the same *corrected* barcode.

4. The file  ``generate_permit_list.json`` that is a JSON file containing information about the run of the command (currently, just the expected orientation).

