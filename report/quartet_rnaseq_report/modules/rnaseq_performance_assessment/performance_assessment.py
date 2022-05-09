#!/usr/bin/env python
""" RnaSeqReport plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging, os
import random
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff

from multiqc.utils import config, report
from multiqc.plots import scatter, table, heatmap
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
        
        # Add to self.css and self.js to be included in template
        self.css = {
            "assets/css/rank.css": os.path.join(
                os.path.dirname(__file__), "assets", "css", "rank.css"
                )
            }
        ########################find log files#################################
        # Plot1: Performance Score
        quality_score = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/quality_score'):
            quality_score = f['f'].splitlines()

        if len(quality_score) != 0:
            ## Now add each section in order
            quality_score_df = self.list2df(quality_score)
            
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
            self.plot_qc_metrics_table(id = 'qc_metrics_summary_table',
                                       qc_metrics_summary_list = table_qc_metrics_dic,
                                       dt_quality_score = quality_score_df)
        else:
            log.debug(
                "No file matched: deg_performance_summary_table - deg_performance_summary.txt"
            )
        print(quality_score_df.loc[quality_score_df['batch'] == 'QC_test', 'rank'].iloc[0])

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
    ##Summary Plot1 and Table1: QC metrics performance summary
    def plot_qc_metrics_table(self,
                              id,
                              qc_metrics_summary_list,
                              dt_quality_score):
        
        # heatmap
        # progress and arrow
        total_len = dt_quality_score.shape[0]
        len_poor = len(dt_quality_score[dt_quality_score["performance"]=='Poor'])
        len_acc = len(dt_quality_score[dt_quality_score["performance"]=='Acceptable'])
        len_out = len(dt_quality_score[dt_quality_score["performance"]=='Outstanding'])
        query_rank = total_len + 1 - int(dt_quality_score.loc[dt_quality_score['batch'] == 'QC_test', 'rank'].iloc[0])
        poor = "%.2f%s" % (len_poor/total_len * 100, '%')
        acceptable = "%.2f%s" % (len_acc/total_len * 100, '%')
        outstanding = "%.2f%s" % (len_out/total_len * 100, '%')
        if query_rank == 1:
            queried = "%.2f%s" % (0, '%')
        elif query_rank == total_len + 1:
            queried = "%.2f%s" % (200, '%')
        else:
            queried = "%.2f%s" % (int(query_rank)/total_len *200, '%')
        
        print({'poor': poor})
        print({'outstanding': outstanding})
        print({'queried': queried})
        
        # ticks number
        snr = dt_quality_score.loc[dt_quality_score['batch'] == 'QC_test', 'SNR'].iloc[0]
        Q0 = dt_quality_score.loc[dt_quality_score['performance'] == 'Poor', 'SNR'].iloc[len_poor -1]
        Q1 = dt_quality_score.loc[dt_quality_score['performance'] == 'Poor', 'SNR'].iloc[0]
        Q2 = dt_quality_score.loc[dt_quality_score['performance'] == 'Acceptable', 'SNR'].iloc[0]
        Q3 = dt_quality_score.loc[dt_quality_score['performance'] == 'Outstanding', 'SNR'].iloc[0]
        print({'Q1': Q1})
        
        # Position of ticks
        tick_Q1 = poor
        tick_Q2 = "%.2f%s" % ((len_poor + len_acc)/total_len * 100, '%')
        print(tick_Q2)
        
        metrics_summary_html = """
        <!-- Arrow -->
        <div class="arrow" style="width: {queried}; margin-top:10px; height: 35px;">
        <svg class="lower-tangle" transform="translate(0, 18)"></svg>
        <span class="lower-label" style="margin-bottom: 25px;"><b> {snr} </b></span>
        </div>
        
        <!-- Progress bar -->
        <div class="progress">
          <div class="progress-bar progress-bar-poor" style="width: {poor}" data-toggle="tooltip" title="" data-original-title="">Poor</div>
          <div class="progress-bar progress-bar-acceptable" style="width: {acceptable}" data-toggle="tooltip" title="" data-original-title="">Acceptable</div>
          <div class="progress-bar progress-bar-outstanding" style="width: {outstanding}" data-toggle="tooltip" title="" data-original-title="">Outstanding</div>
        </div>
        
        <!-- Scale interval -->
        <span style="float:left; left:0%; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q0}</span>
        <span style="float:left; left:{tick_Q1}; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q1}</span>
        <span style="float:left; left:{tick_Q2}; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q2}</span>
        <span style="float:left; left:99%; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q3}</span>
        <br>
        """.format(acceptable=acceptable, outstanding=outstanding, queried=queried, snr=snr, poor=poor, tick_Q1=tick_Q1, tick_Q2=tick_Q2, Q0=Q0, Q1=Q1, Q2=Q2, Q3=Q3)

        # table
        headers = OrderedDict()
        headers['category'] = {
        'title': 'Category',
        'description': 'ategory',
        'scale': False,
        'format': '{0:.0f}'
        }
        headers['value'] = {
        'title': 'Value',
        'description': 'Value',
        'scale': False,
        'format': '{0:.2f}'
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
            'col1_header': 'Quality Metric (s)',
            'no_beeswarm': True,
            'sortRows': False
            }
        
        metrics_table_html = table.plot(qc_metrics_summary_list, headers, table_config)

        # Add a report section with the scatter plot
        self.add_section(
            name='QC metric(s) summary',
            anchor=id + '_anchor',
            description="""
            The submitted data to be tested can be divided into 3 levels based on the SNR by comparing with historical batches: <span style="color: #b80d0d;font-weight:bold">Poor</span>, <span style="color: #70c402;font-weight:bold">Acceptable</span>, <span style="color: #0f9115;font-weight:bold">Outstanding</span>.<br>
            * _Poor_ - the bottom 20%.
            * _Acceptable_ - between bottom 20% and top 80%.
            * _Outstanding_ - the top 20%.
            """,
            plot = metrics_summary_html + '\n' +  metrics_table_html
            )

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
