# Targeted Probe Design Pipeline

## Setup:
Prior to running this pipeline:
* Install all required applications and modules, as listed below.


* Modify the [toml](https://docs.python.org/3.6/library/sqlite3.html) configuration file (probe_design.config.toml) to use on your system.
  * folder path for genome bin fasta files for each cluster
  * folder path for Prokka annotation prediction files (.ffn)
  * Already have an existing blastdb to use?
    * Specify 'use_blastdb' in the configuration file and the pipeline will not create
        new ones for each cluster.


### Application requirements:
* catch/design.py - "A package for designing compact and comprehensive capture probe sets."
  * [https://github.com/broadinstitute/catch/](https://github.com/broadinstitute/catch/)
* ncbi-blast+ for makeblastdb, blastn (used 2.8.1+)
  * [download *blast+* from ncbi.nlm.nih.gov](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
* `clize` required only if pipeline option modification is needed in a separate file.
  * [clize info (ReadTheDocs)](https://clize.readthedocs.io)
  * to install it using pip:  `pip3 install clize`


### Python3:
* [Python](https://www.python.org) &gt;= 3.6
## License
This app is licensed under the terms of the [MIT license](./LICENSE).
