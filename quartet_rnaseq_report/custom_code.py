#!/usr/bin/env python
""" rnaseq-report plugin functions

We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.rnaseq_report_version = get_distribution(
    "quartet_rnaseq_report").version


# Add default config options for the things that are used in MultiQC_NGI
def rnaseq_report_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsedself.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
        return None

    log.info("Running Example MultiQC Plugin v{}".format(
        config.rnaseq_report_version))

    # Add to the main MultiQC config object.
    # User config files have already been loaded at this point
    #   so we check whether the value is already set. This is to avoid
    #   clobbering values that have been customised by users.

    ### Module-rnaseq_data_generation_information
    if 'rnaseq_data_generation_information/information' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_data_generation_information/information': {
                    'fn_re': '^information.json$'
                }
            })

    ### Module-rnaseq_performance_assessment
    if 'rnaseq_performance_assessment/performance_score' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/performance_score': {
                    'fn_re': '^performance_score.txt$'
                }
            })

    if 'rnaseq_performance_assessment/deg_performance_summary' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/deg_performance_summary': {
                    'fn_re': '^deg_performance_summary.txt$'
                }
            })

    if 'rnaseq_performance_assessment/ref_degs_performance_compared' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/ref_degs_performance_compared':
                {
                    'fn_re': '^ref_degs_performance_compared.txt$'
                }
            })

    if 'rnaseq_performance_assessment/rel_performance_summary' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/rel_performance_summary': {
                    'fn_re': '^rel_performance_summary.txt$'
                }
            })

    if 'rnaseq_performance_assessment/ref_rel_exp_per_compared' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/ref_rel_exp_per_compared': {
                    'fn_re': '^ref_rel_exp_per_compared.txt$'
                }
            })

    if 'rnaseq_performance_assessment/studydesign_performance_summary' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/studydesign_performance_summary':
                {
                    'fn_re': '^studydesign_performance_summary.txt$'
                }
            })

    if 'rnaseq_performance_assessment/studydesign_snr' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/studydesign_snr': {
                    'fn_re': '^studydesign_snr.txt$'
                }
            })

    if 'rnaseq_performance_assessment/d5_correlation' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_performance_assessment/d5_correlation': {
                    'fn_re': '^d5_correlation.txt$'
                }
            })
    ### Module-rnaseq_raw_qc
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
    ### Module-post_alignment_qc_modules
    if 'rnaseq_post_alignment_qc/bamqc/genome_results' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/bamqc/genome_results': {
                    'fn': 'genome_results.txt'
                }
            })

    if 'rnaseq_post_alignment_qc/bamqc/coverage' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/bamqc/coverage': {
                    'fn': 'coverage_histogram.txt'
                }
            })

    if 'rnaseq_post_alignment_qc/bamqc/insert_size' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/bamqc/insert_size': {
                    'fn': 'insert_size_histogram.txt'
                }
            })

    if 'rnaseq_post_alignment_qc/bamqc/gc_dist' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/bamqc/gc_dist': {
                    'fn': 'mapped_reads_gc-content_distribution.txt'
                }
            })

    if 'rnaseq_post_alignment_qc/rnaseq/rnaseq_results' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/rnaseq/rnaseq_results': {
                    'fn': 'rnaseq_qc_results.txt'
                }
            })

    if 'rnaseq_post_alignment_qc/rnaseq/coverage' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_post_alignment_qc/rnaseq/coverage': {
                    'fn': 'coverage_profile_along_genes_*'
                }
            })

    ### Module-rnaseq_quantification_qc
    if 'rnaseq_quantification_qc/ref_one_group_performance_summary' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_quantification_qc/ref_one_group_performance_summary': {
                    'fn_re': '^ref_one_group_performance_summary.txt$'
                }
            })

    if 'rnaseq_quantification_qc/ref_detected_gene_performance_compared' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_quantification_qc/ref_detected_gene_performance_compared':
                {
                    'fn_re': '^ref_detected_gene_performance_compared.txt$'
                }
            })

    if 'rnaseq_quantification_qc/ref_two_group_performance_summary' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_quantification_qc/ref_two_group_performance_summary': {
                    'fn_re': '^ref_two_group_performance_summary.txt$'
                }
            })

    if 'rnaseq_quantification_qc/detected_gene_num' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_quantification_qc/detected_gene_num': {
                    'fn_re': '^detected_gene_num.txt$'
                }
            })

    if 'rnaseq_quantification_qc/sd_one_group_performance_summary' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_quantification_qc/sd_one_group_performance_summary': {
                    'fn_re': '^sd_one_group_performance_summary.txt$'
                }
            })

    if 'rnaseq_quantification_qc/sd_cv_mean_ratio' not in config.sp:
        config.update_dict(
            config.sp, {
                'rnaseq_quantification_qc/sd_cv_mean_ratio': {
                    'fn_re': '^sd_one_group_cv.txt$'
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
        'rnaseq_raw_qc', 'rnaseq_post_alignment_qc',
        'rnaseq_quantification_qc', 'rnaseq_supplementary'
    ]
    config.exclude_modules = ['fastqc', 'fastq_screen', 'qualimap']

    config.log_filesize_limit = 2000000000