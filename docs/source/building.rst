Building
========

Alevin-fry is built and tested with the latest (major & minor version) stable Rust (currently 1.45). Building should be as easy as:

.. code:: bash

    $ cargo build --release

subsequent commands below will assume that the executable is in your path. Temporarily, this can be done (in bash-like shells) using:

.. code:: bash

    $ export PATH=`pwd`/target/release/:$PATH

