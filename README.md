# ImmuneBuilder Nextflow Pipeline

<!-- badges: start -->
[![docker-build](https://github.com/juyeongkim/ImmuneBuilder.nf/actions/workflows/docker-build.yml/badge.svg)](https://github.com/juyeongkim/ImmuneBuilder.nf/actions/workflows/docker-build.yml)
<!-- badges: end -->

A Nextflow pipeline based on [ImmuneBuilder](https://github.com/brennanaba/ImmuneBuilder).

## Running Nextflow

### Install Nextflow

Follow this to install Nextflow: https://www.nextflow.io/docs/latest/getstarted.html

### Pull Nextflow pipeline

```sh
nextflow pull juyeongkim/immunebuilder
```

Downloaded pipeline are stored in the folder `$HOME/.nextflow/assets`.

### Run pipeline

#### A) Locally

```sh
cd /where/you/want/to/store/logs/and/intermediate/files
nextflow pull juyeongkim/immunebuilder
nextflow run juyeongkim/immunebuilder -r main --input /your/input/dir --output /your/output/dir
```


#### B) Cluster

Alternatively, you can run the pipeline on HPC with slurm. First, load the environment modules if they are available. If not, please follow the Apptainer and Nextflow documentation to install them first.

```sh
module load Apptainer
module load Nextflow
cd /where/you/want/to/store/logs/and/intermediate/files
nextflow pull juyeongkim/immunebuilder
nextflow run juyeongkim/immunebuilder -r main -profile cluster --input /your/input/dir --output /your/output/dir
```


## Using Docker image

### Get Docker image

```sh
# Pull docker image from GitHub
docker pull ghcr.io/juyeongkim/immunebuilder:latest

# Or build docker image
git clone https://github.com/juyeongkim/immunebuilder.git
cd juyeongkim/immunebuilder.nf
docker build -t ghcr.io/juyeongkim/immunebuilder:latest .
```

### Initiate a Docker container

```sh
cd /WHERE/YOUR/FASTA/FILES/ARE
docker run --volume $PWD:/data -it ghcr.io/juyeongkim/immunebuilder:latest
```

### Run ABodyBuilder2 in Docker container

```sh
cd /data
# your fasta file should have this format: https://github.com/brennanaba/ImmuneBuilder/tree/main#fasta-formatting
ABodyBuilder2 -v --fasta_file yourAntibody.fasta
cat ABodyBuilder2_output.pdb # default outpute file name
```
