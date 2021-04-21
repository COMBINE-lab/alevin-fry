Building
========

Alevin-fry is built and tested with the latest (major & minor) stable version of `Rust <https://www.rust-lang.org/>`_. While it will likely compile fine with older versions of Rust, this is not a guarantee and is not a support priority.  Unlike with C++, Rust has a frequent and stable release cadence, is designed to be installed and updated from user space, and is easy to keep up to date with `rustup <https://rustup.rs/>`_. Thanks to cargo, building should be as easy as:

.. code:: bash

    $ cargo build --release

subsequent commands below will assume that the executable is in your path. Temporarily, this can be done (in bash-like shells) using:

.. code:: bash

    $ export PATH=`pwd`/target/release/:$PATH

