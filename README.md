# Quartet RSeQC Report

The Quartet Project provides publicly accessible multi-omics reference materials and practical tools to enhance the reproducibility and reliability of multi-omics results. Well-characterized multi-omics reference materials and quality control metrics pertinent to precision medicine study purposes can be used to measure and mitigate technical variation, enabling more accurate cross-batch and cross-omics data integration in increasingly large-scale and longitudinal studies such as the International Human Phenome Project. More details on [Quartet Data Portal](https://chinese-quartet.org)

Quartet RSeQC Report is a quality assessment tool for RNA-seq data. It contains two subcommands: `workflow` and `report`. The workflow command takes raw reads (in FASTQ format), produces a set of qc result files from them. and you can use `report` command to report the results finally.

![Workflow](https://docs.chinese-quartet.org/assets/images/rna-workflow.jpeg)

## How to run Quartet RSeQC Report

Assuming your data files are in the /your-data directory.

### 0. Prepare a set of subdirectories

```
mkdir -p /your-data/fastq_screen /your-data/hisat2 /your-data/gtf /your-data/results /your-data/raw-data /your-data/report
```

### 1. Download the dependency files

1. Download reference genomes for `fastq_screen`
wget xxx -O /your-data/fastq_screen

2. Download GTF file
wget xxx -O /your-data/gtf/gencode.v36.annotation.gtf

3. Download Hisat2 index files
wget xxx -O /your-data/hisat2

### 2. Place your data files into `/your-data/raw-data` directory

### 3. Generate qc result files by workflow command

```
docker run -d -v /your-data:/data -it ghcr.io/chinese-quartet/quartet-rseqc-report:latest workflow -i /data/hisat2/GRCh38.d1.vd1.fa.1.ht2 -g /data/gtf/gencode.v36.annotation.gtf -s /data/fastq_screen/fastq_screen.conf --output-dir /data/results --r1 /data/raw-data/example_R1.fq.gz --r2 /data/raw-data/example_R2.fq.gz
```

### 4. Report the results

```
docker run -d -v /your-data:/data -it ghcr.io/chinese-quartet/quartet-rseqc-report:latest report -d /data/results -m /data/metadata.csv --output-dir /data/report
```


## Build from source code

### Prerequisite

- Bash
- Python3 >= 3.7
- pip3
- Java
- R >= 3.6.3

```
conda create -c conda-forge -c bioconda -n quartet-rseqc-report python=3.9 openjdk=8.0.312 r-base=3.6.3
```

### Installation

```
# Activate conda environment
conda activate quartet-rseqc-report

# Clone the repo
git clone https://github.com/chinese-quartet/quartet-rseqc-report

cd quartet-rseqc-report

# Build the environment and compile the quartet-rseqc-report
make all
```

### Usage

```bash
source .env/bin/activate
java -jar target/uberjar/quartet-rseqc-report-*-standalone.jar -h
```

## Plugin Mode

### Prerequisite

Please access [Quartet Service](https://github.com/chinese-quartet/quartet-service) for more details

### Installation

```bash
copm-cli install -n quartet-rseqc-report -V v0.2.2 -d plugins
```

## Examples

...

## Contributions

Quartet QC Report for RNA-Seq Data developed by [Jun Shang](https://github.com/stead99) <[exp2qcdt](./exp2qcdt) & [report](./report)>

## License

Copyright © 2021

This program and the accompanying materials are made available under the
terms of the Eclipse Public License 2.0 which is available at
http://www.eclipse.org/legal/epl-2.0.

This Source Code may also be made available under the following Secondary
Licenses when the conditions for such availability set forth in the Eclipse
Public License, v. 2.0 are satisfied: GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or (at your
option) any later version, with the GNU Classpath Exception which is available
at https://www.gnu.org/software/classpath/license.html.
