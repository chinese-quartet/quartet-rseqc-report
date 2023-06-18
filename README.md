# Quartet RSeQC Report

## How to run Quartet RSeQC Report

see more details on [QDP Docs](https://docs.chinese-quartet.org/data_pipelines/transcriptomics/intro/).


## Build docker image

```bash
bash build-docker.sh

# After build docker image, you can run the docker image
docker run -it --rm quartet-rseqc-report:<tag_name> --help
```

## Build from source code

### Prerequisite

- Bash
- Python3 >= 3.7
- pip3
- Java
- R >= 3.6.3

```
# For Linux
conda create -c conda-forge -c bioconda -c anaconda -n quartet-rseqc-report python=3.9 openjdk=8.0.312 r-base=3.6.3 blas lapack cxx-compiler conda-pack gfortran_linux-64

# For Mac
conda create -c conda-forge -c bioconda -c anaconda -n quartet-rseqc-report python=3.9 openjdk=8.0.312 r-base=3.6.3 blas lapack cxx-compiler conda-pack
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
copm-cli install -n quartet-rseqc-report -V v0.2.4 -d plugins
```

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
