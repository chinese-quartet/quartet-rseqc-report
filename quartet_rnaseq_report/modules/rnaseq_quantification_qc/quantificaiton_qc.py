#!/usr/bin/env python
""" RnaSeqReport plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff

from multiqc import config
from multiqc.plots import scatter, table
from multiqc.modules.base_module import BaseMultiqcModule
from quartet_rnaseq_report.modules.plotly import plot as plotly_plot

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Quantification QC',
            target="quantification_qc",
            anchor='quantification_qc',
            href='https://github.com/clinico-omics/rnaseq-report',
            info=
            " is an report module to show the performance of quartet samples.")

        ### define find log files-------------------------------------------
        # Ref1-table: one group performance summary table based on reference
        ref_one_group_performance_summary = []
        for f in self.find_log_files(
                'rnaseq_quantification_qc/ref_one_group_performance_summary'):
            ref_one_group_performance_summary = f['f'].splitlines()

        if len(ref_one_group_performance_summary) != 0:
            ## Now add each section in order
            ref_one_group_performance_summary_df = self.list2df(
                ref_one_group_performance_summary)
            ref_one_group_performance_summary_df = ref_one_group_performance_summary_df.set_index(
                ['Sample'])
            ref_one_group_performance_summary_dict = ref_one_group_performance_summary_df.to_dict(
                orient='index')
            self.plot_one_group_performance_summary(
                'one_group_performance_summary_table',
                ref_one_group_performance_summary_dict)
        else:
            log.debug(
                "No file matched: one_group_performance_summary_table - plot_one_group_performance_summary.txt"
            )
        # Ref1-plot1: one group performance: section-1 detected gene number performance
        ref_detected_gene_performance_compared = []
        for f in self.find_log_files(
                'rnaseq_quantification_qc/ref_detected_gene_performance_compared'
        ):
            ref_detected_gene_performance_compared = f['f'].splitlines()

        if len(ref_detected_gene_performance_compared) != 0:
            ## Now add each section in order
            ref_detected_gene_performance_compared_df = self.list2df(
                ref_detected_gene_performance_compared)
            self.plot_ref_detected_gene(
                'ref_detected_gene_plot',
                ref_detected_gene_performance_compared_df)
        else:
            log.debug(
                "No file matched: plot_detected_gene_num_bar - detected_gene_num.txt"
            )

        # Ref2-table: two group performance summary table based on reference
        ref_two_group_performance_summary = []
        for f in self.find_log_files(
                'rnaseq_quantification_qc/ref_two_group_performance_summary'):
            ref_two_group_performance_summary = f['f'].splitlines()

        if len(ref_two_group_performance_summary) != 0:
            ## Now add each section in order
            ref_two_group_performance_summary_df = self.list2df(
                ref_two_group_performance_summary)
            ref_two_group_performance_summary_df = ref_two_group_performance_summary_df.set_index(
                ['compared_group'])
            ref_two_group_performance_summary_dict = ref_two_group_performance_summary_df.to_dict(
                orient='index')
            self.plot_two_group_performance_summary(
                'two_group_performance_summary_table',
                ref_two_group_performance_summary_dict)
        else:
            log.debug(
                "No file matched: two_group_performance_summary_table - plot_two_group_performance_summary.txt"
            )

        # Ref2-plot1: two group performance: section-1 relative expression
        ref_rel_exp_per_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/ref_rel_exp_per_compared'):
            ref_rel_exp_per_data = f['f'].splitlines()

        if len(ref_rel_exp_per_data) != 0:
            ## Now add each section in order
            ref_rel_exp_per_data_df = self.list2df(ref_rel_exp_per_data)
            self.plot_ref_rel_exp_per_point(
                'two_group_performance_relative_expression',
                ref_rel_exp_per_data_df)
        else:
            log.debug(
                "No file matched: plot_ref_rel_exp_per_point - ref_rel_exp_per_compared.txt"
            )

        # Ref2-plot2: two group performance: section-1 degs
        ref_degs_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/ref_degs_performance_compared'):
            ref_degs_data = f['f'].splitlines()

        if len(ref_degs_data) != 0:
            ## Now add each section in order
            ref_degs_data_df = self.list2df(ref_degs_data)
            self.plot_ref_degs_point('two_group_performance_degs',
                                     ref_degs_data_df)
        else:
            log.debug(
                "No file matched: ref_degs_point - ref_degs_performance_compared.txt"
            )
        # Study Design Table: summary of detected gene, detected gene JI, CV and correlation value
        sd_one_group_performance_summary = []
        for f in self.find_log_files(
                'rnaseq_quantification_qc/sd_one_group_performance_summary'):
            sd_one_group_performance_summary = f['f'].splitlines()

        if len(sd_one_group_performance_summary) != 0:
            ## Now add each section in order
            sd_one_group_performance_summary_df = self.list2df(
                sd_one_group_performance_summary)
            sd_one_group_performance_summary_df = sd_one_group_performance_summary_df.set_index(
                ['sample'])
            sd_one_group_performance_summary_dict = sd_one_group_performance_summary_df.to_dict(
                orient='index')
            self.plot_sd_one_group_performance_summary(
                'study_design_one_group_performance_summary_table',
                sd_one_group_performance_summary_dict)
        else:
            log.debug(
                "No file matched: study_design_one_group_performance_summary_table - plot_sd_one_group_performance_summary.txt"
            )

        # Study Design plot1:one group performance: section-1 detected gene number in three replicate
        detected_gene_num = []
        for f in self.find_log_files(
                'rnaseq_quantification_qc/detected_gene_num'):
            detected_gene_num = f['f'].splitlines()

        if len(detected_gene_num) != 0:
            ## Now add each section in order
            detected_gene_num_df = self.list2df(detected_gene_num)
            self.plot_detected_gene_num_bar('det_gene_num_plot',
                                            detected_gene_num_df)
        else:
            log.debug(
                "No file matched: plot_detected_gene_num_bar - detected_gene_num.txt"
            )

        # Study Design plot2:one group performance: section-2 CV:mean
        sd_cv_mean_ratio = []
        for f in self.find_log_files(
                'rnaseq_quantification_qc/sd_cv_mean_ratio'):
            sd_cv_mean_ratio = f['f'].splitlines()

        if len(sd_cv_mean_ratio) != 0:
            ## Now add each section in order
            sd_cv_mean_ratio_df = self.list2df(sd_cv_mean_ratio)
            self.plot_cv_mean_point('cv_mean_plot', sd_cv_mean_ratio_df)
        else:
            log.debug(
                "No file matched: plot_cv_mean_point - sd_one_group_cv.txt")

        # Study Design plot3:one group performance: section-3 D5 correlation
        d5_cor_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/d5_correlation'):
            d5_cor_data = f['f'].splitlines()

        if len(d5_cor_data) != 0:
            ## Now add each section in order
            d5_cor_data_df = self.list2df(d5_cor_data)
            self.plot_d5_cor_point('sd_d5_cor_plot', d5_cor_data_df)
        else:
            log.debug("No file matched: snr_pca_point - studydesign_snr.txt")

        # Study Design plot1:more group performance: section-1 SNR
        snr_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/studydesign_snr'):
            snr_data = f['f'].splitlines()

        if len(snr_data) != 0:
            ## Now add each section in order
            snr_data_df = self.list2df(snr_data)
            self.plot_snr_pca_point('sd_mgroup_snr', snr_data_df)
        else:
            log.debug("No file matched: sd_mgroup_snr - studydesign_snr.txt")

    def list2df(self, data):
        """Convert string list to dataframe"""
        cols = data[0].split('\t')
        array = []
        for line in data[1:]:
            array.append(line.split('\t'))
        df = pd.DataFrame(np.array(array), columns=cols)
        return df

    ### Add sections-----------------------------------------------------
    # One group summary table: summary detected gene between reference and this batch with a table
    def plot_one_group_performance_summary(
            self,
            id,
            ref_one_group_performance_summary_dict,
            title="Detected Gene Summary",
            section_name=None,
            description=None,
            helptext=None):
        headers = OrderedDict()
        headers['Sensitivity'] = {
            'title': 'Sensitivity',
            'description':
            'Sensitivity is the proportion of "true" detected genes from reference dataset which can be correctly detected by the test set.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['Reference_Sensitivity_min'] = {
            'title': 'Reference (Min_of_Sensitivity)',
            'description':
            'Sensitivity is the proportion of "true" detected genes from reference dataset which can be correctly detected by the test set.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['Reference_Sensitivity_max'] = {
            'title': 'Reference (Max_of_Sensitivity)',
            'description':
            'Sensitivity is the proportion of "true" detected genes from reference dataset which can be correctly detected by the test set.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['Specificity'] = {
            'title': 'Specificity',
            'description':
            'Specificity is the proportion of "true" non-detected genes from reference dataset which can be correctly not detected by the test set.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['Reference_Specificity_min'] = {
            'title': 'Reference (Max_of_Specificity)',
            'description':
            'Specificity is the proportion of "true" non-detected genes from reference dataset which can be correctly not detected by the test set.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['Reference_Specificity_max'] = {
            'title': 'Reference (Max_of_Specificity)',
            'description':
            'Specificity is the proportion of "true" non-detected genes from reference dataset which can be correctly not detected by the test set.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        self.add_section(
            name='A group summarize  based on reference data',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNAs',
            helptext=helptext if helptext else '''
                This longer description explains what exactly the numbers mean
                and supports markdown formatting. This means that we can do _this_:

                * Something important
                * Something else important
                * Best of all - some `code`

                Doesn't matter if this is copied from documentation - makes it
                easier for people to find quickly.
                ''',
            plot=table.plot(ref_one_group_performance_summary_dict, headers))

    # One group performance compare: comare detected gene between reference and this batch
    def plot_ref_detected_gene(
            self,
            id,
            ref_degs_data_df,
            title="Reference-OneGroup-Detected Gene Performance",
            section_name=None,
            description=None,
            helptext=None):
        fig = px.scatter(ref_degs_data_df,
                         x="Sensitivity_tier1",
                         y="Specificity_tier1",
                         color="Group",
                         marginal_x="box",
                         marginal_y="box",
                         template="simple_white")

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='Reference-OneGroup-Detected Gene Performance',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNAs',
            helptext=helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot=html)

    # Two group summary table: relative expression and degs performance based on reference
    def plot_two_group_performance_summary(
            self,
            id,
            ref_two_group_performance_summary_dict,
            title="Two Group Performance Summary",
            section_name=None,
            description=None,
            helptext=None):
        headers = OrderedDict()
        headers['consistent'] = {
            'title': 'Consistency ratio of relative expression',
            'description':
            'Proportion of genes that falls into reference range (mean ± 2 fold SD) in relative ratio (log2FC).',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['consistent_min'] = {
            'title': 'Reference (Min_of_Consistent)',
            'description': 'Min reference of consistent',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['consistent_max'] = {
            'title': 'Reference (Max_of_Consistent)',
            'description': 'Max reference of consistent',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['corr'] = {
            'title': 'Correlation of relative log2FC',
            'description':
            'Pearson correlation between mean value of reference relative ratio and test site.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['corr_min'] = {
            'title': 'Reference (Min_of_Correlation)',
            'description': 'Min reference of Correlation',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['corr_max'] = {
            'title': 'Reference (Max_of_Correlation)',
            'description': 'Max reference of Correlation',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['Sensitivity'] = {
            'title': 'Sensitivity of DEGs',
            'description':
            'Sensitivity is the proportion of "true" DEGs from reference dataset which can be correctly identified as DEG by the test set.',
            'max': 1,
            'min': 0,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['sensitivity_min'] = {
            'title': 'Reference (Min_of_Sensitivity)',
            'description': 'Min reference of Sensitivity.',
            'max': 1,
            'min': 0,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['sensitivity_max'] = {
            'title': 'Reference (Max_of_Sensitivity)',
            'description': 'Max reference of Sensitivity.',
            'max': 1,
            'min': 0,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }
        headers['Specificity'] = {
            'title': 'Specificity of DEGs',
            'description':
            'Specificity is the proportion of "true" not DEGs from reference dataset which can be can be correctly identified as non-DEG by the test set.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['specificity_min'] = {
            'title': 'Reference (Min_of_Specificity)',
            'description': 'Min reference of Specificity.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['specificity_max'] = {
            'title': 'Reference (Max_of_Specificity)',
            'description': 'Max reference of Specificity.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        self.add_section(
            name='Two group summarize  based on reference data',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNAs',
            helptext=helptext if helptext else '''
                This longer description explains what exactly the numbers mean
                and supports markdown formatting. This means that we can do _this_:

                * Something important
                * Something else important
                * Best of all - some `code`

                Doesn't matter if this is copied from documentation - makes it
                easier for people to find quickly.
                ''',
            plot=table.plot(ref_two_group_performance_summary_dict, headers))

    # Two group performance compare: relative expression performance based on reference
    def plot_ref_rel_exp_per_point(
            self,
            id,
            ref_rel_exp_per_data_df,
            title="Reference-TwoGroup-Relative Expression Performance",
            section_name=None,
            description=None,
            helptext=None):
        fig = px.scatter(ref_rel_exp_per_data_df,
                         x="corr",
                         y="consistent",
                         color="Group",
                         marginal_x="box",
                         marginal_y="box",
                         template="simple_white")

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='Reference-TwoGroup-Relative Expression Performance',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNAs',
            helptext=helptext if helptext else '''
                This longer description explains what exactly the numbers mean
                and supports markdown formatting. This means that we can do _this_:

                * Something important
                * Something else important
                * Best of all - some `code`

                Doesn't matter if this is copied from documentation - makes it
                easier for people to find quickly.
                ''',
            plot=html)

    # Two group performance compare: degs performance based on reference
    def plot_ref_degs_point(self,
                            id,
                            ref_degs_data_df,
                            title="Reference-TwoGroup-DEG Performance",
                            section_name=None,
                            description=None,
                            helptext=None):
        fig = px.scatter(ref_degs_data_df,
                         x="Sensitivity",
                         y="Specificity",
                         color="Group",
                         marginal_x="box",
                         marginal_y="box",
                         template="simple_white")

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='Reference-TwoGroup-DEG Performance',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNAs',
            helptext=helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot=html)

    # Study Design table: summary of detected gene, detected gene JI, CV and correlation value
    def plot_sd_one_group_performance_summary(
            self,
            id,
            sd_one_group_performance_summary_dict,
            title="StudyDesign-OneGroup-Summary",
            section_name=None,
            description=None,
            helptext=None):
        headers = OrderedDict()
        headers['N'] = {
            'title': 'Number of detected genes',
            'description':
            'This metric is used to estimate the detection abundance of one sample.',
            'max': 25000,
            'min': 15000,
            'suffix': '',
            'scale': 'RdYlGn',
            'format': '{:,.2f}'
        }

        headers['ref_n_min'] = {
            'title': 'Reference (Min_of_DetectedGene)',
            'description': 'Min Reference of N.',
            'max': 25000,
            'min': 15000,
            'suffix': '',
            'scale': 'RdYlGn',
            'format': '{:,.2f}'
        }

        headers['ref_n_max'] = {
            'title': 'Reference (Max_of_DetectedGene)',
            'description': 'Max Reference of N.',
            'max': 25000,
            'min': 15000,
            'suffix': '',
            'scale': 'RdYlGn',
            'format': '{:,.2f}'
        }

        headers['JacardIndex'] = {
            'title': 'Detection Jaccard index',
            'description':
            'Detection JI is the ratio of number of the genes detected in both replicates than the number of the genes detected in either of the replicates. This metric is used to estimate the repeatability of one sample detected gene from different replicates.',
            'max': 1,
            'min': 0.6,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['ref_ji_min'] = {
            'title': 'Reference (Min_of_JI)',
            'description': 'Min reference of JI',
            'max': 1,
            'min': 0.6,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['ref_ji_max'] = {
            'title': 'Reference (Max_of_JI)',
            'description': 'Max reference of JI',
            'max': 1,
            'min': 0.6,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['CV'] = {
            'title': 'Correlation of relative log2FC',
            'description':
            'CV is calculated based on the normalized expression levels in all 3 replicates of one sample for each genes. This metric is used to estimate the repeatability of one sample expression level from different replicates.',
            'max': 1000000000,
            'min': 5,
            'suffix': '',
            'scale': 'OrRd',
            'format': '{:,.2f}'
        }

        headers['ref_cv_min'] = {
            'title': 'Reference (Min_of_CV)',
            'description': 'Min reference of CV',
            'max': 1000000000,
            'min': 5,
            'suffix': '',
            'scale': 'OrRd',
            'format': '{:,.2f}'
        }

        headers['ref_cv_max'] = {
            'title': 'Reference (Max_of_CV)',
            'description': 'Max reference of CV',
            'max': 1000000000,
            'min': 5,
            'suffix': '',
            'scale': 'OrRd',
            'format': '{:,.2f}'
        }

        headers['CTR'] = {
            'title': 'Correlation of technical replicates',
            'description':
            'CTR is calculated based on the correlation of one sample expression level from different replicates.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['ref_ctr_min'] = {
            'title': 'Reference (Min_of_CTR)',
            'description': 'Min reference of CTR.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        headers['ref_ctr_max'] = {
            'title': 'Reference (Max_of_CTR)',
            'description': 'Max reference of CTR.',
            'max': 1,
            'min': 0.8,
            'suffix': '',
            'scale': 'Blue',
            'format': '{:,.2f}'
        }

        self.add_section(
            name='StudyDesign-OneGroup-Summary',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNAs',
            helptext=helptext if helptext else '''
                This longer description explains what exactly the numbers mean
                and supports markdown formatting. This means that we can do _this_:

                * Something important
                * Something else important
                * Best of all - some `code`

                Doesn't matter if this is copied from documentation - makes it
                easier for people to find quickly.
                ''',
            plot=table.plot(sd_one_group_performance_summary_dict, headers))

    # Study Design plot1:one group performance: section-1 detected gene number in three replicate
    def plot_detected_gene_num_bar(self,
                                   id,
                                   detected_gene_num_df,
                                   title="Reference DEGs Performance Compared",
                                   section_name=None,
                                   description=None,
                                   helptext=None):
        fig = px.bar(detected_gene_num_df,
                     x="sample",
                     y="N",
                     color="n_rep",
                     title="Detected Gene Number",
                     template="simple_white")

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='Relative Expression Performance Based On Reference',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNAs',
            helptext=helptext if helptext else '''
                This longer description explains what exactly the numbers mean
                and supports markdown formatting. This means that we can do _this_:

                * Something important
                * Something else important
                * Best of all - some `code`

                Doesn't matter if this is copied from documentation - makes it
                easier for people to find quickly.
                ''',
            plot=html)

    # Study Design plot2:one group performance: section-2 cv correlation
    def plot_cv_mean_point(self,
                           id,
                           sd_cv_mean_ratio_df,
                           title="Study Design-One Group Performance-CV",
                           section_name=None,
                           description=None,
                           helptext=None):
        fig = px.scatter(sd_cv_mean_ratio_df,
                         x="mean",
                         y="CV",
                         template="simple_white")
        fig.update_layout(title='The distribution of intra-batch CV ',
                          xaxis_title='Mean FPKM',
                          yaxis_title='CV',
                          font=dict(family="Courier New, monospace",
                                    size=18,
                                    color="black"))

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='Study Design-One Group Performance-CV',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNAs',
            helptext=helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot=html)

    # Study Design plot3:one group performance: section-3 d5 correlation
    def plot_d5_cor_point(
            self,
            id,
            d5_cor_data_df,
            title="Study Design-One Group Performance-Correlation of expression",
            section_name=None,
            description=None,
            helptext=None):
        fig = px.scatter(d5_cor_data_df,
                         x="D5_1",
                         y="D5_2",
                         template="simple_white")
        fig.update_layout(title='CTR',
                          xaxis_title='D5_1',
                          yaxis_title='D5_2',
                          font=dict(family="Courier New, monospace",
                                    size=18,
                                    color="black"))

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='Study Design-One Group Performance-Correlation of expression',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNAs',
            helptext=helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot=html)

    # Study Design plot1:more group performance: section-1 SNR
    def plot_snr_pca_point(self,
                           id,
                           snr_data_df,
                           title="Study Design-More Group Performance-SNR",
                           section_name=None,
                           description=None,
                           helptext=None):
        fig = px.scatter(snr_data_df,
                         x="PC1",
                         y="PC2",
                         color="sample",
                         color_discrete_map={
                             "D5": "#4CC3D9",
                             "D6": "#7BC8A4",
                             "F7": "#FFC65D",
                             "M8": "#F16745"
                         })
        SNR_value = 'SNR (' + snr_data_df.iloc[1].at['SNR'] + ')'
        PC1_value = 'PC1 (' + snr_data_df.iloc[1].at['PC1_ratio'] + '%)'
        PC2_value = 'PC2 (' + snr_data_df.iloc[1].at['PC2_ratio'] + '%)'
        fig.update_layout(title=SNR_value,
                          xaxis_title=PC1_value,
                          yaxis_title=PC2_value,
                          legend_title='Sample',
                          font=dict(family="Arial, sans-serif",
                                    size=18,
                                    color="black"),
                          template="simple_white")

        fig.update_xaxes(showline=True,
                         linewidth=1,
                         linecolor='black',
                         mirror=True)

        fig.update_yaxes(showline=True,
                         linewidth=1,
                         linecolor='black',
                         mirror=True)

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='Study Design-More Group Performance-SNR',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNAs',
            helptext=helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot=html)
