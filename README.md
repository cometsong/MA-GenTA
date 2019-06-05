# Targeted Probe Design Pipeline

## Setup:
Prior to running this pipeline:
* Install all required applications and modules, as listed below.

* Modify the [toml](https://github.com/toml-lang/toml/) configuration file (targeted_probe_config.toml) or make a new one to use on your system.
  * Set folder path for genome bin fasta files for each cluster.  These paths can be relative to the folder you start the pipeline in, or absolute from root filesystem.
  * Set folder path for Prokka annotation prediction files (.ffn) ... *OR*
  * Already have an existing blastdb to use?
    * Specify 'use_blastdb' in your configuration file and the pipeline will not create
        new ones for each cluster.


### RunningPipeline
After setting up folders and configuration file, call the pipeline:  
  > `python3 <path/to/script/>targeted_probe_design.py`

or for modified configuration options:  
  > `python3 <path/to/script/>targeted_probe_design.py --config-file <filename.toml>`  


### Application requirements:
- **catch/design.py** &gt;= 1.2.0
  * "A package for designing compact and comprehensive capture probe sets."
  * [github.com/broadinstitute/catch](https://github.com/broadinstitute/catch/)
- **ncbi-blast+** for makeblastdb, blastn (used 2.8.1+)
  * [ftp.ncbi.nlm.nih.gov](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
- **sqlite3** database client
  * [sqlite.org](https://sqlite.org/)
- **Python3**
  * [&gt;= v3.6](https://www.python.org/downloads/)
  * modules:
    - **logbook**: for logging to screen and file ( `pip3 install logbook` )
    - **tomlkit**: used for all config options    ( `pip3 install tomlkit` )
    - **clize** required *only* if pipeline option modification is needed in a separate file.
      * [clize info (ReadTheDocs)](https://clize.readthedocs.io)
      * to install: `pip3 install clize`
      * Once installed, add `--config-file <filename.toml>`
        to your command line for run-specific options

## License
This app is licensed under the terms of the [MIT license](./LICENSE).
