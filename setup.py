#!/usr/bin/env python
"""
MultiReport for Quartet RNAseq Report.
"""

from setuptools import setup, find_packages

version = '0.1.1'

setup(
    name='quartet_rnaseq_report',
    version=version,
    author='Jun Shang',
    author_email='shangjunv@163.com',
    description="MultiReport for Quartet RNA-Seq Pipeline.",
    long_description=__doc__,
    keywords='bioinformatics',
    url='https://github.com/clinico-omics/quartet_rnaseq_report',
    download_url=
    'https://github.com/clinico-omics/quartet_rnaseq_report/releases',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['multiqc==1.9', 'plotly==4.9.0', 'pandas==1.1.0'],
    entry_points={
        'multiqc.modules.v1': [
            'rnaseq_data_generation_information = rnaseq_report.modules.rnaseq_data_generation_information:MultiqcModule',
            'rnaseq_performance_assessment = rnaseq_report.modules.rnaseq_performance_assessment:MultiqcModule',
            'rnaseq_raw_qc = rnaseq_report.modules.rnaseq_raw_qc:MultiqcModule',
            'rnaseq_post_alignment_qc = rnaseq_report.modules.rnaseq_post_alignment_qc:MultiqcModule',
            'rnaseq_quantification_qc = rnaseq_report.modules.rnaseq_quantification_qc:MultiqcModule',
            'rnaseq_supplementary = rnaseq_report.modules.rnaseq_supplementary:MultiqcModule'
            # 'rnaseq_qc = rnaseq_report.modules.rnaseq_qc:MultiqcModule'
        ],
        'multiqc.cli_options.v1':
        ['disable_plugin = quartet-rnaseq-report.cli:disable_plugin'],
        'multiqc.hooks.v1': [
            'execution_start = rnaseq_report.custom_code:rnaseq_report_execution_start'
        ],
        'multiqc.templates.v1': ['default = rnaseq_report.templates.default']
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
