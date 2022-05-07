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
        # Plot1: Performance Score
        quality_score = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/quality_score'):
            quality_score = f['f'].splitlines()

        if len(quality_score) != 0:
            ## Now add each section in order
            quality_score_df = self.list2df(quality_score)
            print(quality_score_df)
            quality_score_df.astype({'quality_score': float})
            test_score = quality_score_df.loc[
                quality_score_df['batch'] ==
                'QC_test']['SNR'].values[0]
            quality_score_list = [[
                float(i)
                for i in quality_score_df['SNR'].values.tolist()
            ]]
            print(test_score)
            print(quality_score_list)
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
            table_qc_metrics_dic = qc_metrics_summary_df.set_index('qc_metrics').T.to_dict()
            print(table_qc_metrics_dic)
            self.plot_qc_metrics_table(id = 'qc_metrics_summary_table',
                                       qc_metrics_summary_list = table_qc_metrics_dic)
        else:
            log.debug(
                "No file matched: deg_performance_summary_table - deg_performance_summary.txt"
            )

        # Plot2: SNR
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


    def list2df(self, data):
        """Convert string list to dataframe"""
        cols = data[0].split('\t')
        array = []
        for line in data[1:]:
            array.append(line.split('\t'))
        df = pd.DataFrame(np.array(array), columns=cols)
        return df

    #######################add section#################################
    ###Plot1: quality score
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
            name='Quality score',
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
                              description=None,
                              helptext=None):
        headers = OrderedDict()
        headers['category'] = {
        'title': 'Category',
        'description': 'ategory',
        'scale': False,
        'format': '{0:.3f}'
        }
        headers['value'] = {
        'title': 'Value',
        'description': 'Value',
        'scale': False,
        'format': '{0:.3f}'
        }
        headers['historical_value'] = {
            'title': 'Historical value (mean ± SD)',
            'description': 'Historical Value (mean ± SD)',
            'scale': False,
            'format': '{0:.2f}'
            }
        headers['rank'] = {
            'title': 'Rank',
            'description': 'Rank',
            'scale': False,
            'format': '{:.0f}'
            }
        table_config = {
            'namespace': 'conclusion_summary',
            'id': id,
            'table_title': '',
            'col1_header': 'Quality Metric(s)',
            'no_beeswarm': True,
            'sortRows': False
            }
        
        metrics_html = table.plot(qc_metrics_summary_list, headers, table_config)

        # Add a report section with the scatter plot
        self.add_section(
            name='QC metrics summary',
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
            plot=metrics_html)

    # Plot2: SNR performance evaluation
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
