Getting Started
===============

There are a (growing) number of different sub-commands: generate-permit-list, collate, and quant. Each of these is invoked as a command passed as the first argument to the alevin-fry executable. For example, to run the generate-permit-list command, one would run:

.. code:: bash

    $ alevin-fry generate-permit-list --help

This should then show the following:

.. code:: bash

    $ alevin-fry generate-permit-list --help
    alevin-fry-generate-permit-list 0.0.1
    Avi Srivastava, Rob Patro
    Generate a permit list of barcodes from a RAD file

    USAGE:
        alevin-fry generate-permit-list [FLAGS] --input <input> --output-dir <output-dir> --expect-cells <expect-cells> --force-cells <force-cells> --valid-bc <valid-bc>

    FLAGS:
        -h, --help             Prints help information
        -k, --knee-distance    attempt to determine the number of barcodes to keep using the knee distance method
        -V, --version          Prints version information

    OPTIONS:
        -e, --expect-cells <expect-cells>    defines the expected number of cells to use in determining the (read, not UMI) based cutoff
        -f, --force-cells <force-cells>      select the top-k most-frequent barcodes, based on read count, as valid (true)
        -i, --input <input>                  input RAD file
        -o, --output-dir <output-dir>        output directory
        -b, --valid-bc <valid-bc>            uses true barcode collected from a provided file
