#!/usr/bin/env python3

import os
import re
import json
import click
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
              type=click.Path(exists=True, file_okay=True),
              help="The fastq file with suffixes of _R1.fastq.gz or _R1.fq.gz.")
@click.option('--r2', required=True,
              type=click.Path(exists=True, file_okay=True),
              help="The fastq file with suffixes of _R2.fastq.gz or _R2.fq.gz.")
@click.option('--hisat2-index', '-i', required=True,
              type=click.Path(exists=True, dir_okay=True),
              help="The index for the reference genome.")
@click.option('--gtf', '-g', required=True,
              type=click.Path(exists=True, file_okay=True),
              help="The gtf file.")
@click.option('--output-dir', required=False,
              type=click.Path(exists=True, dir_okay=True),
              help="The output directory.")
@click.option('--fastq-screen-conf', '-s', required=True,
              type=click.Path(exists=True, file_okay=True),
              help="The config file for fastq-screen, the reference genomes must be located in the same directory with config file.")
def workflow(r1, r2, hisat2_index, fastq_screen_conf, gtf, output_dir):
    if not re.match(r'_R1.(fastq|fq).gz', r1):
        raise Exception(
            "The fastq file must be with suffixes of _R1.fastq.gz or _R1.fq.gz")

    if not re.match(r'_R2.(fastq|fq).gz', r2):
        raise Exception(
            "The fastq file must be with suffixes of _R2.fastq.gz or _R2.fq.gz")

    wdl_dir = '/venv/workflow'
    if not os.path.exists(wdl_dir):
        print("Cannot find the workflow, please contact the administrator.")

    def zip_path_by_ext_program(input_path, output_path):
        cmd = ['zip', '-r', '-q', output_path, input_path]
        print('ZIP: Working Directory %s, CMD: %s' % (os.getcwd(), cmd))
        proc = Popen(cmd, stdin=PIPE)
        proc.communicate()

    zip_path_by_ext_program('/venv/workflow/tasks', "/venv/workflow/tasks.zip")

    data_dict = {
        "rseqc.read1": r1,
        "rseqc.read2": r2,
        "rseqc.idx": hisat2_index,
        "rseqc.fastq_screen_conf": fastq_screen_conf,
        "rseqc.gtf": gtf
    }

    inputs_file = "/venv/workflow/inputs"
    inputs = read_json(inputs_file)
    inputs.update(data_dict)

    write_json(inputs, inputs_file)

    def call_cromwell(inputs_fpath, workflow_fpath, workflow_root, tasks_path):
        cmd = ['cromwell', 'run', workflow_fpath, "-i", inputs_fpath, "-p", tasks_path, "--workflow-root", workflow_root]
        print('Run workflow and output results to %s.' % workflow_root)
        proc = Popen(cmd, stdin=PIPE)
        proc.communicate()

    inputs_fpath = "/venv/workflow/inputs"
    workflow_fpath = "/venv/workflow/workflow.wdl"
    workflow_root = output_dir
    tasks_path = "/venv/workflow/tasks.zip"
    call_cromwell(inputs_fpath, workflow_fpath, workflow_root, tasks_path)