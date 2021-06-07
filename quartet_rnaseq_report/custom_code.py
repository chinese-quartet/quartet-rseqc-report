#!/usr/bin/env python
""" rnaseq-report plugin functions

We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.quartet_rnaseq_report_version = get_distribution("quartet_rnaseq_report").version


# Add default config options for the things that are used in MultiQC_NGI
def quartet_rnaseq_report_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsedself.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
        return None

    log.info("Running Example MultiQC Plugin v{}".format(
        config.quartet_rnaseq_report_version))

    # Add to the main MultiQC config object.
    # User config files have already been loaded at this point
    #   so we check whether the value is already set. This is to avoid
    #   clobbering values that have been customised by users.

    # Module-rnaseq_data_generation_information
    if 'rnaseq_data_generation_information/information' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_data_generation_information/information': {
                    'fn_re': '^general-info.json$'
                }
            })

    # Module-rnaseq_performance_assessment
    if 'rnaseq_performance_assessment/quality_score' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/quality_score': {
                    'fn_re': '^quality_score.txt$'
                }
            })

    if 'rnaseq_performance_assessment/performance_of_absolute_exp' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/performance_of_absolute_exp': {
                    'fn_re': '^performance_of_absolute_exp.txt$'
                }
            })

    if 'rnaseq_performance_assessment/performance_of_relative_exp' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/performance_of_relative_exp': {
                    'fn_re': '^performance_of_relative_exp.txt$'
                }
            })

    if 'rnaseq_performance_assessment/pca_with_snr' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/pca_with_snr': {
                    'fn_re': '^pca_with_snr.txt$'
                }
            })

    if 'rnaseq_performance_assessment/relative_exp_correlation' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/relative_exp_correlation': {
                    'fn_re': '^relative_exp_correlation.txt$'
                }
            })

    if 'rnaseq_performance_assessment/absolute_exp_correlation' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/absolute_exp_correlation': {
                    'fn_re': '^absolute_exp_correlation.txt$'
                }
            })

    if 'rnaseq_performance_assessment/qc_metrics_summary' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/qc_metrics_summary': {
                    'fn_re': '^qc_metrics_summary.txt$'
                }
            })
    # Module-rnaseq_raw_qc
    if 'rnaseq_raw_qc/zip' not in config.sp:
        config.update_dict(config.sp,
                           {'rnaseq_raw_qc/zip': {
                               'fn': '*_fastqc.zip'
                           }})

    if 'rnaseq_raw_qc/data' not in config.sp:
        config.update_dict(config.sp,
                           {'rnaseq_raw_qc/data': {
                               'fn': 'fastqc_data.txt'
                           }})

    if 'rnaseq_raw_qc/fastq_screen' not in config.sp:
        config.update_dict(
            config.sp, {'rnaseq_raw_qc/fastq_screen': {
                'fn': '*_screen.txt'
            }})
    # Module-post_alignment_qc_modules
    if 'rnaseq_post_alignment_qc/bam_qc/genome_results' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/bam_qc/genome_results': {
                    'fn': 'genome_results.txt'
                }
            })

    if 'rnaseq_post_alignment_qc/bam_qc/coverage' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/bam_qc/coverage': {
                    'fn': 'coverage_histogram.txt'
                }
            })

    if 'rnaseq_post_alignment_qc/bam_qc/insert_size' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/bam_qc/insert_size': {
                    'fn': 'insert_size_histogram.txt'
                }
            })

    if 'rnaseq_post_alignment_qc/bam_qc/gc_dist' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/bam_qc/gc_dist': {
                    'fn': 'mapped_reads_gc-content_distribution.txt'
                }
            })

    if 'rnaseq_post_alignment_qc/rnaseq_qc/rnaseq_qc_results' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/rnaseq_qc/rnaseq_qc_results': {
                    'fn': 'rnaseq_qc_results.txt'
                }
            })

    if 'rnaseq_post_alignment_qc/rnaseq_qc/coverage' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/rnaseq_qc/coverage': {
                    'fn': 'coverage_profile_along_genes_*'
                }
            })

    # # Some additional filename cleaning
    # config.fn_clean_exts.extend([
    #     '.my_tool_extension',
    #     '.removeMetoo'
    # ])

    # # Ignore some files generated by the custom pipeline
    # config.fn_ignore_paths.extend([
    #     '*/my_awesome_pipeline/fake_news/*',
    #     '*/my_awesome_pipeline/red_herrings/*',
    #     '*/my_awesome_pipeline/noisy_data/*',
    #     '*/my_awesome_pipeline/rubbish/*'
    # ])

    config.module_order = [
        'rnaseq_data_generation_information', 'rnaseq_performance_assessment',
        'rnaseq_raw_qc', 'rnaseq_post_alignment_qc', 'rnaseq_supplementary'
    ]
    config.exclude_modules = ['fastqc', 'fastq_screen', 'qualimap']

    config.log_filesize_limit = 2000000000
