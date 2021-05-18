# single-cell
Benchmarking variant calling in single cells

### Running experiments

All experiments are reproducable using a Snakemake workflow. First clone the repository:

```shell
$ git clone https://github.com/luntergroup/single-cell && cd single-cell
```

You will need to install Conda, if not already:

```shell
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh # follow instructions, answer 'yes' where asked
$ source ~/.bashrc # assuming you installed conda into your home directory
$ conda update conda
```

Install Snakemake and general dependencies with conda:

```shell
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda create --name single-cell snakemake pysam python-wget openpyxl
$ conda activate single-cell
```

Each set of experiments is specified in a YAML config file in the `config` directory. For example, to run the fibroblast analysis:

```shell
$ snakemake --configfile config/dong.yaml --use-conda -j 100 --cluster "qsub -cwd -V -j y -P mygroup.prj -q long.qf -pe shmem {threads}"
```
