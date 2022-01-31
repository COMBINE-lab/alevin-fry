Other resources for alevin-fry
==============================

In addition to the current documentation page, there are numerous other resources to help you learn more about alevin-fry, how to process data using
this program, and how to further process the output of alevin-fry in downstream analysis.

Tutorials
---------

A collection of tutorials describing how to process different types of data with `alevin-fry` and describing different features of `alevin-fry` is 
available `here <https://combine-lab.github.io/alevin-fry-tutorials/#blog>`_.

FAQ
---

We hope to make use of GitHub discussions to answer frequently asked questions, and to discuss other issues relevant to the development and use
of `alevin-fry`.  You can visit the GitHub discussion page for `alevin-fry here <https://github.com/COMBINE-lab/alevin-fry/discussions>`_.  
GitHub discussions are also a good place to raise large-scale feature requests to see if they make sense in the context of `alevin-fry`.  For 
small-scale feature requests, or to report bugs or unexpected behavior you encounter when processing data with `alevin-fry`, please make use 
of our `GitHub issues page <https://github.com/COMBINE-lab/alevin-fry/issues>`_.

Quality Control
---------------

Support for `alevin-fry` in the `alevinQC <https://github.com/csoneson/alevinQC>`_ package is imminent.

Easy loading of USA-mode data
-----------------------------

The `fishpond <https://mikelove.github.io/fishpond/>`_ package contains many methods for making the ingestion of quantification results generated 
by `salmon <https://github.com/COMBINE-lab/salmon>`_ and `alevin-fry` into R easy.  In particular, you can find documentation on the 
`loadFry function here <https://mikelove.github.io/fishpond/reference/loadFry.html>`_.  This makes it easy to import USA-mode quantification 
results into a `SingleCellExperiment <https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html>`_ object, and to properly 
extract or combine the spliced, unspliced, and ambiguous count components.