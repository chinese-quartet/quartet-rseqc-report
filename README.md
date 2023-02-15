# quartet-rseqc-report

Visualizes Quality Control(QC) results for Quartet Project.

## Standalone Mode

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

## Run with Docker

```
docker run -d -v /root/rseqc-raw.py:/venv/bin/rseqc.py -v /root/cromwell-local.conf:/venv/cromwell-local.conf -v /mnt/home_ssd/home/yangjingcheng/test_quartet_rseqc_report:/data -v /mnt/home_ssd:/mnt/home_ssd -v /root/workflow-raw:/venv/workflow -v /mnt/pgx_src_data_pool_4:/mnt/pgx_src_data_pool_4 -it ghcr.io/chinese-quartet/quartet-rseqc-report:7eaaa39-7eaaa395 workflow -i /mnt/pgx_src_data_pool_4/reference/human/GRCh38/hisat2/GRCh38.d1.vd1.fa.1.ht2 -g /mnt/pgx_src_data_pool_4/reference/human/GRCh38/annotation_files/gencode.v36.annotation.gtf -s /mnt/pgx_src_data_pool_4/reference/human/GRCh38/fastq_screen/fastq_screen.conf --output-dir /mnt/home_ssd/home/yangjingcheng/test_quartet_rseqc_report --r1 /mnt/pgx_src_data_pool_4/fuscc_lc_1000/rnaseq/clean/2568LC_R1.fq.gz --r2 /mnt/pgx_src_data_pool_4/fuscc_lc_1000/rnaseq/clean/2568LC_R2.fq.gz
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
