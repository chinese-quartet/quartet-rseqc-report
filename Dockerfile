###################
# STAGE 1: builder
###################

# Build currently doesn't work on > Java 11 (i18n utils are busted) so build on 8 until we fix this
FROM debian:stable-slim as builder

WORKDIR /app/source

ENV PATH="$PATH:/opt/conda/bin:/opt/conda/envs/venv/bin"
ENV FC_LANG en-US
ENV LC_CTYPE en_US.UTF-8

# bash:    various shell scripts
# wget:    installing lein
# git:     ./bin/version
# make:    backend building
# gettext: translations
RUN apt-get update && apt-get install -y coreutils bash git wget make gettext

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_22.11.1-1-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p /opt/conda
RUN /opt/conda/bin/conda install -c conda-forge -c bioconda -c anaconda mamba blas lapack cxx-compiler conda-pack gfortran_linux-64

# Note: cromwell==83 must not deleted.
RUN /opt/conda/bin/mamba create -n venv -c bioconda -c conda-forge -y python=3.9 r-renv r-base=3.6.3 hisat2==2.2.1 samtools bioconductor-ballgown bioconductor-genefilter qualimap==2.2.2d fastq-screen==0.15.2 fastqc==0.11.9 fastp==0.23.2 stringtie==2.2.1 cromwell==83

# Customized softewares
ADD ./bin/quartet-rseqc-report /opt/conda/envs/venv/bin/quartet-rseqc-report
ADD ./bin/rseqc.py /opt/conda/envs/venv/bin/rseqc.py

# Install report locally instead of remote to install the latest version.
ADD report /report
RUN /opt/conda/envs/venv/bin/pip install /report

ADD ./resources/bin/exp2qcdt.sh /opt/conda/envs/venv/bin/exp2qcdt.sh
ADD ./resources/renv /opt/conda/envs/venv/renv
ADD ./resources/renv.lock /opt/conda/envs/venv/renv.lock
# Install exp2qcdt locally instead of remote to install the latest version.
# ADD ./build/Rprofile /opt/conda/envs/venv/etc/Rprofile
# RUN /opt/conda/envs/venv/bin/Rscript /opt/conda/envs/venv/etc/Rprofile

# Disable cache to install all packages into the conda environment.
COPY exp2qcdt /exp2qcdt
RUN /opt/conda/envs/venv/bin/Rscript -e 'renv::activate("/opt/conda/envs/venv");renv::restore();renv::install("/exp2qcdt")'

# For app render.
RUN /opt/conda/envs/venv/bin/pip install git+https://github.com/yjcyxky/biominer-app-util.git

# Build quartet-rseqc-report.jar
# Add the rest of the source
ADD . .
# Fetch all submodule
RUN git submodule update --init --recursive

# lein: backend dependencies and building
ADD ./bin/lein /usr/local/bin/lein
RUN chmod 744 /usr/local/bin/lein

# install dependencies before adding the rest of the source to maximize caching
# backend dependencies
ADD project.clj .
RUN lein deps

# build the app
RUN lein uberjar

# Pack the conda environment
RUN conda-pack -n venv -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

RUN /venv/bin/conda-unpack

# ###################
# # STAGE 2: runner
# ###################

FROM debian:stable-slim as runner

ENV PATH="$PATH:/venv/bin"
ENV PYTHONDONTWRITEBYTECODE=1
ENV FC_LANG en-US
ENV LC_CTYPE en_US.UTF-8

LABEL org.opencontainers.image.source = "https://github.com/chinese-quartet/quartet-rseqc-report"

RUN apt-get update && apt-get install -y coreutils bash git wget

WORKDIR /data

# Customized
COPY --from=builder /venv /venv
COPY --from=builder /app/source/target/uberjar/quartet-rseqc-report*.jar /quartet-rseqc-report.jar
COPY --from=builder /app/source/resources/Rprofile /venv/etc/Rprofile

# Copy all the installed R packages from the builder.
COPY --from=builder /root/.local/share/renv /root/.local/share/renv

# Initialize Rprofile when running the container. The src/quartet_rseqc_report/cli.clj will check the Rprofile to find the path of the renv environment.
RUN sed -i 's/<plugin_env_path>/\/venv/g' /venv/etc/Rprofile

## Workflow - WDL files
COPY --from=builder /app/source/workflow /venv/workflow

## Config file for cromwell instance
COPY --from=builder /app/source/build/cromwell-local.conf /venv/cromwell-local.conf

## Make count work properly.
RUN ln -s /venv/bin/prepDE.py /venv/bin/count

## Add ballgown wrapper
COPY ./build/ballgown /venv/bin/ballgown
RUN chmod a+x /venv/bin/ballgown

# Run it
ENTRYPOINT ["rseqc.py"]