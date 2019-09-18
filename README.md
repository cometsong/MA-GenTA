# Targeted Probe Design Pipeline

----------

## Overview:
> Design targeted assay probes using `prokka` annotated prediction files,
> genome bins files, `catch`, and `blast` databases.

This script serves as a pipeline to design unique probes for a genome (or
group of genomes) and customize the selection of each probe based on
user-defined specifications.

These probes can be used for example in the detection and relative
quantification of bacteria in microbiome samples using high throughput
sequencing.

----------

## Base Setup:

Prior to running this pipeline, install all required applications and modules,
as listed below, then modifying all configuration settings to match your
system.

#### Application requirements:
For each of the apps listed, either ensure they are in the `$PATH` variable on
your system, or use the `[APPS]` section of the config file to set their paths.

- **catch/design.py** &gt;= 1.2.0  ([github.com/broadinstitute/catch][catch])
  * *"A package for designing compact and comprehensive capture probe sets."*
- **ncbi-blast+** for makeblastdb, blastn &gt;=2.8.1  ([ftp.ncbi.nlm.nih.gov][blast])
- **sqlite3** database client &gt;=3.7.17  ([sqlite.org][])
- **Python3** [&gt;= v3.6][py3]
  * modules:
    - **logbook**: for logging to screen and file ( `pip3 install logbook` )
    - **tomlkit**: used for all config options    ( `pip3 install tomlkit` )
    - **clize**: required _only_ if using new pipeline config file path on the
        command line. ( `pip3 install clize` )
      * To use: add option `--config-file <awesome-config-file.toml>`  
        to your command line for run-specific configuration

Using a `virtualenv` is _highly_ recommended. :)

----------

## Running Pipeline:

#### Configure Your Variables
* Modify the default [toml][] configuration file ([targeted_probe_config.toml][config])  
  _or_ copy to a new one to use on your system.
  - `[paths]`:  
    These paths can be relative to the folder you start the pipeline in, or
        absolute in the root filesystem.
    * working_dir: folder to hold resulting files
    * genome_bins: folder path for fasta files
    * prokka_dir: for Prokka annotation prediction files (.ffn) from which
        individual blastdb files will be created.
    * use_blastdb: see _issue_ below for info on this
  - `[general]`:
    * final_probe_amount: you want to produce (default 20)
    * final_probe_random: choose randomly among resulting probes? (default true)
    * prokka_prediction_suffix: suffix on prediction files (default '.ffn')
    * genome_bins_suffix: suffix of bins files (default '.fasta')
  - `[qc_percent]`:
    * min_percent: default 45
    * max_percent: default 65
  - `[catch]`:
    * probe_length: how long resulting probes? (default 40)
    * probe_stride: how  many base pairs between probes? (default 20)
    * reuse_existing_probe_files: reuse preexisting results, e.g. a previous
        pipeline run failed (default false)
  - `[blastn]`:
    * evalue: for cutoff of blast resulting records (default '0.001')
    * num_alignments: Integer >1. (blastn default: 250)
    * num_threads: how many cpus use? (default 2)
    * fields: add extra fields to the default set [qseqid, sseqid, pident, length, qseq]
  - `[filters]`:
    * musicc_list: set of strings to match for results to _keep_
    * trna_list: set of strings to match for results to _skip_
  - `[APPS]`:
    * See above for more info on this section.

#### Run the Script
- Using the default configuration file, call the pipeline:
  > `python3 targeted_probe_design.py`

- or after setting options in a modified configuration file:
  > `python3 targeted_probe_design.py --config-file awesome-config-file.toml`

#### Results
The resulting files from each run of this pipeline will include:
- fasta file containing sequences of filtered matching probes 
- tsv file with output of _coverage analysis_ from running `catch`
- modified version of the prokka prediction file (header includes ffn name with
    with spaces replaced by underscores)
- sqlite database with all matching probe info and sequences in a table, and
    a view of the probes filtered according to the _config_ file settings
- a log file

The probes fasta file is the one you want. 
All others can be reviewed, kept for records, gzip'd, or trashed as you see fit.

----------

### Issues? FAQ!
  * _Already have an existing blastdb to use:_
    * Specify `use_blastdb = 'path/to/db'` in `[paths]` section of your 
        configuration file and that will be used for all pipeline runs
        instead of creating new ones for each `.ffn` file.

  * _Dealing with plethora of genomes?_
    * If processing the entire set of files in a single pipeline run causes
        issues such as too looooong (or too much *walltime* on your computing
        cluster), then splitting it into chunks is one method.
        * Create _boilerplate_ config file; setup the `[paths]` section of the
          config file so it can be _serialized_ by incrementing numbers:
          ```
          [paths]
               genome_bins=/path/to/gbins/01  # example base numbered dirname
               # or if splitting prokka annotation files:
               prokka_dir=/path/to/ffns/01
          ```
          saved with a filename like `n01.toml`.
        * Make serially-named directories; move or link (using `ln`) each set
          genome bins (or annotated files...) into each dir.
        * Use `seq` and `sed` to modify and create config files for each dir,
          this example creating sequential files `n02.toml` to `n15.toml`:
          ```
          line="genome_bins=";
          cur=01; max=15;
          for D in $(seq -w $cur $((max-1)) ); do
            sed -r 's/(.*)(${line}.*)([0-9]+)(.*)/echo "\1\2$((\3+1))\4"/ge' \
            n${D}.toml \
            > n$(printf "%0${#cur}d" $((D+1))).toml ;
          done
          ``` 
          Setting the `cur` to the number in the _boilerplate_ filename, and the
          `max` to the highest directory number. Any leading _zero_'s on the
          `cur` will be retained in the example above.

----------

## License
This app is licensed under the terms of the MIT [license][].


[LINKS]:'reference-list'
[config]:./targeted_probe_config.toml
[toml]:https://github.com/toml-lang/toml/blob/master/README.md
[catch]:https://github.com/broadinstitute/catch/blob/master/README.md#catch-----
[blast]:https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/
[sqlite.org]:https://sqlite.org/
[py3]:https://www.python.org/downloads/
[license]:./LICENSE
