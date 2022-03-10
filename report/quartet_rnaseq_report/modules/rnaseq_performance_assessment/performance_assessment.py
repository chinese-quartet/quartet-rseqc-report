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
from multiqc.plots import scatter
from multiqc.modules.base_module import BaseMultiqcModule
from quartet_rnaseq_report.modules.plotly import plot as plotly_plot

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Performance Assessment',
            target="performance_assessment",
            anchor='performance_assessment',
            href='https://github.com/clinico-omics/quartet-rnaseq-report',
            info=
            " is an report module to show the performance of quartet samples.")
        ########################find log files#################################
        # Plot0: Performance Score
        quality_score = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/quality_score'):
            quality_score = f['f'].splitlines()

        if len(quality_score) != 0:
            ## Now add each section in order
            quality_score_df = self.list2df(quality_score)
            quality_score_df.astype({'quality_score': float})
            test_score = quality_score_df.loc[
                quality_score_df['batch'] ==
                'QC_test']['quality_score'].values[0]
            quality_score_list = [[
                float(i)
                for i in quality_score_df['quality_score'].values.tolist()
            ]]
            self.plot_quality_score('plot_quality_score', quality_score_list,
                                    test_score)
        else:
            log.debug("No file matched: quality_score - quality_score.txt")

        # Summary Table1: ALL QC metrics summary
        qc_metrics_summary = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/qc_metrics_summary'):
            qc_metrics_summary = f['f'].splitlines()

        if len(qc_metrics_summary) != 0:
            ## Now add each section in order
            qc_metrics_summary_df = self.list2df(qc_metrics_summary)
            qc_metrics_summary_list = [
                qc_metrics_summary_df.columns.values.tolist()
            ] + qc_metrics_summary_df.values.tolist()
            self.plot_qc_metrics_table('qc_metrics_summary_table',
                                       qc_metrics_summary_list)
        else:
            log.debug(
                "No file matched: deg_performance_summary_table - deg_performance_summary.txt"
            )

        # Plot1: performance of absolute expression
        dt_abs_exp_evaluate = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/performance_of_absolute_exp'):
            dt_abs_exp_evaluate = f['f'].splitlines()

        if len(dt_abs_exp_evaluate) != 0:
            ## Now add each section in order
            abs_exp_evaluate_df = self.list2df(dt_abs_exp_evaluate)
            self.plot_abs_exp_evaluate_point('abs_exp_evaluate_plot',
                                             abs_exp_evaluate_df)
        else:
            log.debug(
                "No file matched: performance_of_absolute - performance_of_absolute_exp.txt"
            )

        # Plot2: performance of relative expression
        dt_rel_exp_evaluate = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/performance_of_relative_exp'):
            dt_rel_exp_evaluate = f['f'].splitlines()

        if len(dt_rel_exp_evaluate) != 0:
            ## Now add each section in order
            rel_exp_evaluate_df = self.list2df(dt_rel_exp_evaluate)
            self.plot_rel_exp_evaluate_point('el_exp_evaluate_plot',
                                             rel_exp_evaluate_df)
        else:
            log.debug(
                "No file matched: performance_of_relative - performance_of_relative_exp.txt"
            )

        # Plot3: SNR
        snr_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/pca_with_snr'):
            snr_data = f['f'].splitlines()

        if len(snr_data) != 0:
            ## Now add each section in order
            snr_data_df = self.list2df(snr_data)
            self.plot_snr_pca_point('snr_pca_plot', snr_data_df)
        else:
            log.debug("No file matched: snr_pca_point - studydesign_snr.txt")

        # Plot4: Relative correlation
        rel_exp_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/relative_exp_correlation'):
            rel_exp_data = f['f'].splitlines()

        if len(rel_exp_data) != 0:
            ## Now add each section in order
            rel_exp_data_df = self.list2df(rel_exp_data)
            self.plot_rel_cor_point('rel_cor_plot', rel_exp_data_df)
        else:
            log.debug(
                "No file matched: rel_cor_plot - relative_exp_correlation.txt")

        # Plot5: Abosolute correlation
        abs_exp_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/absolute_exp_correlation'):
            abs_exp_data = f['f'].splitlines()

        if len(abs_exp_data) != 0:
            ## Now add each section in order
            abs_exp_data_df = self.list2df(abs_exp_data)
            self.plot_abs_cor_point('abs_cor_plot', abs_exp_data_df)
        else:
            log.debug(
                "No file matched: abs_cor_plot - absolute_exp_correlation.txt")

    def list2df(self, data):
        """Convert string list to dataframe"""
        cols = data[0].split('\t')
        array = []
        for line in data[1:]:
            array.append(line.split('\t'))
        df = pd.DataFrame(np.array(array), columns=cols)
        return df

    #######################add section#################################
    ###Plot0: quality score
    def plot_quality_score(self,
                           id,
                           quality_score_list,
                           test_score,
                           title=None,
                           section_name=None,
                           description=None,
                           helptext=None):
        fig = px.imshow(quality_score_list,
                        x=quality_score_list[0],
                        y=['score'],
                        template="simple_white")
        fig.update_traces(dict(showscale=False,
                               coloraxis=None,
                               colorscale='RdYlGn'),
                          selector={'type': 'heatmap'})
        fig.update_layout(showlegend=False,
                          annotations=[
                              dict(x=test_score,
                                   y=-2.5,
                                   xref="x",
                                   yref="y",
                                   text=str(test_score)[0:4],
                                   textangle=0,
                                   showarrow=True,
                                   font=dict(family="Arial, sans-serif",
                                             size=45,
                                             color="black"),
                                   align="center",
                                   arrowhead=1,
                                   arrowsize=4,
                                   arrowwidth=3,
                                   arrowside="end",
                                   arrowcolor="grey",
                                   ax=0,
                                   ay=-60,
                                   yshift=-145)
                          ])
        fig.update_xaxes(ticks="outside",
                         tickwidth=2,
                         tickcolor='black',
                         ticklen=10,
                         showline=True,
                         linewidth=2,
                         linecolor='black',
                         tickfont=dict(family='Arial, sans-serif',
                                       color='black',
                                       size=20))
        fig.update_yaxes(showline=False, showticklabels=False, showgrid=False)

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='Quality Score',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using reference RNA',
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

    ##Table1: QC metrics performance summary
    def plot_qc_metrics_table(self,
                              id,
                              qc_metrics_summary_list,
                              title="QC Metrics Summary",
                              section_name=None,
                              description=None,
                              helptext=None):
        fig = ff.create_table(qc_metrics_summary_list, height_constant=60)

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='QC Metrics Summary',
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

    # Plot1: Abosolute expression performance evaluation
    def plot_abs_exp_evaluate_point(
            self,
            id,
            abs_exp_evaluate_df,
            title="Performance evaluation on the intra-batch level",
            section_name=None,
            description=None,
            helptext=None):
        fig = px.scatter(abs_exp_evaluate_df,
                         x="SNR",
                         y="LIR",
                         color="protocol",
                         color_discrete_map={
                             "P": "#2f5c85",
                             "R": "#7ba1c7",
                             "QC": "red"
                         },
                         marginal_x="box",
                         marginal_y="box",
                         template="simple_white")

        fig.update_layout(xaxis_title='SNR',
                          yaxis_title='Absolute correlation',
                          font=dict(family="Arial, sans-serif",
                                    size=18,
                                    color="black"),
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
            name='Performance evaluation on the intra-batch level',
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

    # Plot2: Relative expression performance evaluation
    def plot_rel_exp_evaluate_point(
            self,
            id,
            rel_exp_evaluate_df,
            title="Performance evaluation in relative expression",
            section_name=None,
            description=None,
            helptext=None):
        fig = px.scatter(rel_exp_evaluate_df,
                         x="corr_ref",
                         y="corr_FC",
                         color="protocol",
                         color_discrete_map={
                             "P": "#2f5c85",
                             "R": "#7ba1c7",
                             "QC": "red"
                         },
                         marginal_x="box",
                         marginal_y="box",
                         template="simple_white")
        fig.update_layout(
            xaxis_title='Reference datasets based on relative correlation',
            yaxis_title='Relative correlation',
            font=dict(family="Arial, sans-serif", size=18, color="black"),
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
            name='Performance evaluation in relative expression',
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

    # Plot3: SNR performance evaluation
    def plot_snr_pca_point(self,
                           id,
                           snr_data_df,
                           title=None,
                           section_name=None,
                           description=None,
                           helptext=None):
        SNR_value = 'SNR = ' + snr_data_df.iloc[1].at[
            'SNR'] + ' (N = ' + snr_data_df.iloc[1].at['gene_num'] + ')'
        PC1_value = 'PC1 (' + snr_data_df.iloc[1].at['PC1_ratio'] + '%)'
        PC2_value = 'PC2 (' + snr_data_df.iloc[1].at['PC2_ratio'] + '%)'
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
        fig.update_traces(marker=dict(size=16))
        fig.update_layout(title=SNR_value,
                          xaxis_title=PC1_value,
                          yaxis_title=PC2_value,
                          legend_title='Sample',
                          font=dict(family="Arial, sans-serif",
                                    size=18,
                                    color="black"),
                          template="simple_white")
        fig.update_layout(title={
            'y': 0.99,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        })
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
            name='Signal-to-Noise Ratio',
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

    # Plot4: Relative Correlation figure
    def plot_rel_cor_point(self,
                           id,
                           rel_exp_data_df,
                           title=None,
                           section_name=None,
                           description=None,
                           helptext=None):
        pccs_rel = str(
            round(
                rel_exp_data_df.iloc[:, 0].astype(float).corr(
                    rel_exp_data_df.iloc[:, 1].astype(float)), 3))
        gene_num = str(len(rel_exp_data_df))
        fig = px.scatter(rel_exp_data_df,
                         x=rel_exp_data_df.iloc[:, 0],
                         y=rel_exp_data_df.iloc[:, 1],
                         color_discrete_sequence=['#7BC8A4'],
                         hover_name='gene_id',
                         template="simple_white")
        fig.update_layout(title='cor = ' + pccs_rel + ' (N = ' + gene_num +
                          ')',
                          xaxis_title=rel_exp_data_df.columns[0],
                          yaxis_title=rel_exp_data_df.columns[1],
                          font=dict(family="Arial, sans-serif",
                                    size=18,
                                    color="black"))
        fig.update_layout(title={
            'y': 0.99,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        })
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
            name='Relative Correlation',
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

    # Plot5: Abosolute Correlation figure

    def plot_abs_cor_point(self,
                           id,
                           abs_exp_data_df,
                           title=None,
                           section_name=None,
                           description=None,
                           helptext=None):
        pccs_abs = str(
            round(
                abs_exp_data_df.iloc[:, 1].astype(float).corr(
                    abs_exp_data_df.iloc[:, 2].astype(float)), 3))
        gene_num = str(len(abs_exp_data_df))
        fig = px.scatter(abs_exp_data_df,
                         x=abs_exp_data_df.iloc[:, 1],
                         y=abs_exp_data_df.iloc[:, 2],
                         hover_name='gene_id',
                         template="simple_white")
        fig.update_layout(title='cor = ' + pccs_abs + ' (N = ' + gene_num +
                          ')',
                          xaxis_title=abs_exp_data_df.columns[1],
                          yaxis_title=abs_exp_data_df.columns[2],
                          font=dict(family="Arial, sans-serif",
                                    size=18,
                                    color="black"))
        fig.update_layout(title={
            'y': 0.99,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        })
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
            name='Abosolute Correlation',
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
