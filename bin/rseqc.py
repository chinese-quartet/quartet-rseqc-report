#!/usr/bin/env python3

import os
import re
import json
import click
# You may need to install https://github.com/yjcyxky/biominer-app-util firstly.
from biominer_app_util.cli import render_app
from subprocess import Popen, PIPE


def read_json(json_file):
    with open(json_file, "r") as f:
        return json.load(f)


def write_json(data, json_file):
    with open(json_file, 'w') as f:
        json.dump(data, f)


@click.group()
def rseqc():
    pass


@rseqc.command(help="Run the pipeline for RNA-Seq data.")
@click.option('--r1', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help="The fastq file with suffixes of _R1.fastq.gz or _R1.fq.gz.")
@click.option('--r2', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help="The fastq file with suffixes of _R2.fastq.gz or _R2.fq.gz.")
@click.option('--hisat2-index', '-i', required=True,
              type=click.Path(exists=True, dir_okay=True, file_okay=False),
              help="The index for the reference genome.")
@click.option('--gtf', '-g', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help="The gtf file.")
@click.option('--output-dir', required=False,
              type=click.Path(exists=True, dir_okay=True, file_okay=False),
              help="The output directory.")
@click.option('--fastq-screen-conf', '-s', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help="The config file for fastq-screen, the reference genomes must be located in the same directory with config file.")
def workflow(r1, r2, hisat2_index, fastq_screen_conf, gtf, output_dir):
    if not re.match(r'.*_R1.(fastq|fq).gz', r1):
        raise Exception(
            "The R1 fastq file must be with suffixes of _R1.fastq.gz or _R1.fq.gz")

    if not re.match(r'.*_R2.(fastq|fq).gz', r2):
        raise Exception(
            "The R2 fastq file must be with suffixes of _R2.fastq.gz or _R2.fq.gz")

    wdl_dir = '/venv/workflow'
    if not os.path.exists(wdl_dir):
        print("Cannot find the workflow, please contact the administrator.")

    for item in ["GRCh38.d1.vd1.fa.1.ht2", "GRCh38.d1.vd1.fa.2.ht2", "GRCh38.d1.vd1.fa.3.ht2",
              "GRCh38.d1.vd1.fa.4.ht2", "GRCh38.d1.vd1.fa.5.ht2", "GRCh38.d1.vd1.fa.6.ht2",
              "GRCh38.d1.vd1.fa.7.ht2", "GRCh38.d1.vd1.fa.8.ht2"]:
        if not os.path.exists(os.path.join(hisat2_index, item)):
            raise Exception("Cannot find %s in %s, you need to download hisat2 index files." % (item, hisat2_index))

    project_name = "rseqc"
    data_dict = {
        "project_name": project_name,
        "read1": r1,
        "read2": r2,
        "idx": hisat2_index,
        "fastq_screen_conf": fastq_screen_conf,
        "gtf": gtf
    }

    output_workflow_dir = os.path.join(os.path.dirname(output_dir), "workflow")
    os.makedirs(output_workflow_dir, exist_ok=True)

    render_app(wdl_dir, output_dir=output_workflow_dir,
               project_name=project_name, sample=data_dict)

    def call_cromwell(inputs_fpath, workflow_fpath, workflow_root, tasks_path):
        # cmd = ['cromwell', 'run', workflow_fpath, "-i", inputs_fpath,
        #        "-p", tasks_path, "--workflow-root", workflow_root]
        cmd = ['java', '-Dconfig.file=/venv/cromwell-local.conf', '-jar', '/venv/share/cromwell/cromwell.jar', 'run', workflow_fpath, "-i", inputs_fpath, "-p", tasks_path, "--workflow-root", workflow_root]
        print('Run workflow and output results to %s.' % workflow_root)
        proc = Popen(cmd, stdin=PIPE)
        proc.communicate()

    inputs_fpath = os.path.join(output_workflow_dir, "inputs")
    workflow_fpath = os.path.join(output_workflow_dir, "workflow.wdl")
    tasks_path = os.path.join(output_workflow_dir, "tasks.zip")
    call_cromwell(inputs_fpath, workflow_fpath, output_dir, tasks_path)


@rseqc.command(help="Run the report for RNA-Seq results.")
@click.option('--result-dir', '-d', required=True,
              type=click.Path(exists=True, file_okay=True),
              help="A directory which contains a series of results from RNA-Seq pipeline.")
@click.option('--metadata-file', '-m', required=True,
              type=click.Path(exists=True, file_okay=True),
              help="A metadata file")
@click.option('--output-dir', '-o', required=True,
              type=click.Path(exists=True, dir_okay=True),
              help="A directory which will store the output report.")
def report(result_dir, metadata_file, output_dir):
    cmd = ['quartet-rseqc-report', '-d', result_dir, "-m", metadata_file,
           "-o", output_dir]
    print('Run quartet-rseqc-report and output the report to %s.' % output_dir)
    proc = Popen(cmd, stdin=PIPE)
    proc.communicate()


if __name__ == '__main__':
    rseqc()
