Installing alevin-fry
=====================

Alevin-fry can be installed using a package manager such as ``conda``, or built from source.

Installing with bioconda
------------------------

Alevin-fry is available for both x86 linux and OSX platforms `using bioconda <https://anaconda.org/bioconda/alevin-fry>`_.

With ``bioconda`` in the appropriate place in your channel list, you should simply be able to install via:

.. code:: bash

    $ conda install alevin-fry

Installing from source
----------------------

If you want to use features or fixes that may only be available in the latest develop branch (or want to build for a different 
architecture), then you have to build from source.  Luckily, ``cargo`` makes that easy; see below.

Alevin-fry is built and tested with the latest (major & minor) stable version of `Rust <https://www.rust-lang.org/>`_. While it will likely compile fine with older versions of Rust, this is not a guarantee and is not a support priority.  Unlike with C++, Rust has a frequent and stable release cadence, is designed to be installed and updated from user space, and is easy to keep up to date with `rustup <https://rustup.rs/>`_. Thanks to cargo, building should be as easy as:

.. code:: bash

    $ cargo build --release

subsequent you will want to place ``alevin-fry`` in your ``PATH``. This can be done (in bash-like shells) using:

.. code:: bash

    $ export PATH=`pwd`/target/release/:$PATH

To ensure that ``alevin-fry`` remains in your path between logins, you should make sure the path to ``target/release/`` shown above is set in the ``PATH`` variable in the appropriate file for your shell (e.g. in ``~/.profile``, ``~/.bashrc`` etc.).
