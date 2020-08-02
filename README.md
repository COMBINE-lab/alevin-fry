# alevin-fry

Absolutely RADic(a)l methods for analyzing single-cell sequencing data (written in Rust)!

## Building

Alevin-fry is built and tested with the latest (major & minor version) stable Rust (currently 1.45).
Building should be as easy as:

```{bash}
$ cargo build --release
```

subsequent commands below will assume that the executable is in your path.  Temporarily, this can 
be done (in bash-like shells) using:

```{bash}
$ export PATH=`pwd`/target/release/:$PATH
```

## Running alevin-fry

There are a (growing) number of different sub-commands: `generate-permit-list`, `collate`, and `quant`.  Each of these
is invoked as a command passed as the first argument to the `alevin-fry` executable.  For example, to run the 
`generate-permit-list` command, one would run:

```{bash}
$ alevin-fry generate-permit-list --help
```

This should then show the following:

```{bash}
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
```

The commands and options are described below.


### generate-permit-list

This command takes as input a RAD file (created by running `alevin` with the `--justAlign` flag), and it determines what cell barcodes should
be associated with "true" cells, which should be corrected to some "true" barcode, and which should simply be ignored / discarded.  This command 
has 3 required arguments; the path to an input RAD file `--input`, the path to an output directory `--output-dir` (which will be created if it doesn't exist),
and then **one** of the following mutually exclusive options (which determines how the "true" barcodes are decided):

  * `--knee-distance` : This flag will use the `distance` method that is used in the `whitelist` command of [UMI-tools](https://github.com/CGATOxford/UMI-tools) to attempt to automatically determine the number of true barcodes.  Briefly, this method first counts the number of _reads_ associated with each barcode, and then sorts the barcodes in descending order by their associated read count.  It then constructs the cumulative distribution function from this sorted list of frequencies.  Finally, it applies an iterative algorithm to attempt to determine the optimal number of barcodes to include by looking for a "knee" or "elbow" in the CDF graph.  The algorithm considers each barcode in the CDF where it's x-coordinate is equal to this barcode's rank divided by the total number of barcodes (i.e. its normalized rank) and the y-coordinate is equal to the (normalized) cumulative frequency achieved at this barcode.  It then computes the distance of this barcode from the line x=y (defined by the start and end of the CDF).  The initial knee is predicted as the point that has the maximum distance from the x=y line.  The algorithm is iterative, because experiments with many low-quality barcodes may predict too many valid barcodes using this method.  Thus, the algorithm is run repeatedly, each time considering a prefix of the CDF from index 0 through the previous knee's index * 5.  Once two subsequent iterations of the algorithm return the same knee point, the algorithm terminates.
  * `--force-cells <ncells>`: This option will count the number of reads associated with each barcode, and sort the barcodes in descending order of frequency.  Then, it will consider the first `<ncells>` barcodes to be valid.  Any barcode that has a number of reads >= to the `<ncells>`-th barcode will be considered part of the permit list, all others will not (but will be considered for _correction_ to this permit list).
  * `--valid-bc <bcfile>`: This option will read the provided file `<bcfile>` and treat it as an explicitly-provided list of true barcodes.  Barcodes appearing in this list will be considered true, and barcodes will be corrected to this list.
  * `--expect-cells <ncells>`: Not currently implemnted.
  
  

