Getting Started
===============

Running the alevin-fry pipeline
-------------------------------

First, we need to generate a RAD file using alevin.  The RAD file is created by mapping the sequencing reads against an index of the reference. We recommend using a `splici <https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/>`_ reference index. The mappings can be generated using either `selective-alignment <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8>`_ or `pseudoalignment <https://www.nature.com/articles/nbt.3519>`_ against the transcriptome (with the ``--rad`` or ``--sketch`` flags, respectively). **Note**, however, that alevin-fry does not currently support RAD files aligned against a decoy-aware index, so that indices used for RAD file generation should be prepared without decoy sequnece. For a chromium v2 set of read files, the command would look like the following:

.. code:: bash

    $ salmon alevin -lISR --chromium -1 <read1_files> -2 <read2_files> -o <alevin_odir> -i <index> -p <num_threads> --sketch

Given the output directory generated above, the next step is to let alevin-fry generate the permit list.  First, we grab the 10x Chromium version 2
permit list (if we had Chromium v3 chemistry, we would use that permit-list instead):

.. code:: bash

    $ wget https://umd.box.com/shared/static/jbs2wszgbj7k4ic2hass9ts6nhqkwq1p -O 10x_v2_permit.txt

Now, we can use this permit list to scan the cell barcodes actually encountered in our reads and determine a set of cells that were likely present in our sample:

.. code:: bash 

    $ alevin-fry generate-permit-list --input <alevin_odir> --expected-ori fw --output-dir <fry_odir> --unfiltered-pl 10x_v2_permit.txt

Next, given the permit list and barcode mapping (which resides in the `<fry_odir>` directory), we collate the original RAD file using the command below.

.. code:: bash 

    $ alevin-fry collate -i <fry_odir> -r <alevin_odir> -t <num_threads>

Finally, we quantify the collated rad file using the `cr-like` resolution strategy using the `quant` command below.

.. code:: bash 

    $ alevin-fry quant -i <fry_odir> -m <tg_map> -t <num_threads> -r cr-like -o <fry_odir> 

Note that with the exception of the `generate-permit-list` command, the other `alevin-fry` commands are designed to scale well with the number of provided threads. Thus, if you have multiple threads to use, then you can provide the appropriate argument to the `-t` option.

Detailed information on the alevin-fry commands
-----------------------------------------------

There are a (growing) number of different sub-commands for ``alevin-fry``.  To learn more about the different commands an their options check the :ref:`commands<alevin-fry commands>` section of the documentation.
