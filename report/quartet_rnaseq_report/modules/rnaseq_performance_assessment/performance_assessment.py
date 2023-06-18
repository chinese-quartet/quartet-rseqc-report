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
            href='https://github.com/chinese-quartet/quartet-rseqc-report',
            info=
            " is an report module to show the performance of quartet samples.")
        
        # Add to self.css and self.js to be included in template
        self.css = {
            "assets/css/rank.css": os.path.join(
                os.path.dirname(__file__), "assets", "css", "rank.css"
                )
            }
        ########################find log files#################################
        ### Quality Score Table
        quality_score = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/quality_score'):
            quality_score = f['f'].splitlines()

        if len(quality_score) != 0:
            ## Now add each section in order
            dt_quality_score = self.list2df(quality_score)
            
        else:
            log.debug("No file matched: quality_score - quality_score.txt")

        # Summary Table1: ALL QC metrics summary
        qc_metrics_summary = []
        for f in self.find_log_files('rnaseq_performance_assessment/qc_metrics_summary'):
            qc_metrics_summary = f['f'].splitlines()
            
        if len(qc_metrics_summary) != 0:
            ## Now add each section in order
            qc_metrics_summary_df = self.list2df(qc_metrics_summary)
            table_qc_metrics_dic = qc_metrics_summary_df.set_index('qc_metrics').T.to_dict()
            self.plot_qc_metrics_table(id = 'qc_metrics_summary_table',
                                       qc_metrics_summary_list = table_qc_metrics_dic,
                                       dt_quality_score = dt_quality_score)
        else:
            log.debug(
                "No file matched: deg_performance_summary_table - deg_performance_summary.txt"
            )
        
        # Plot2: SNR and RC scatter plot
        self.plot_snr_rc_point('plot_snr_rc_point', dt_quality_score)

        # Plot3: SNR
        list_snr = []
        for f in self.find_log_files('rnaseq_performance_assessment/pca_with_snr'):
            list_snr = f['f'].splitlines()

        if len(list_snr) != 0:
            ## Now add each section in order
            df_snr = self.list2df(list_snr)
            self.plot_snr_pca_point('snr_pca_plot', df_snr)
        else:
            log.debug("No file matched: snr_pca_point - studydesign_snr.txt")
        
        # Plot4: RC correlation
        list_rc = []
        for f in self.find_log_files('rnaseq_performance_assessment/logfc_cor_ref_test'):
            list_rc = f['f'].splitlines()
            
        if len(list_rc) != 0:
            ## Now add each section in order
            df_rc = self.list2df(list_rc)
            self.plot_rc_cor_point('rc_cor_plot', df_rc)
        else:
            log.debug(
                "No file matched: performance_of_relative - logfc_cor_ref_test.txt"
            )

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
        len_bad = len(dt_quality_score[dt_quality_score["performance"]=='Bad'])
        len_fair = len(dt_quality_score[dt_quality_score["performance"]=='Fair'])
        len_good = len(dt_quality_score[dt_quality_score["performance"]=='Good'])
        len_great = len(dt_quality_score[dt_quality_score["performance"]=='Great'])
        query_rank = total_len + 1 - int(dt_quality_score.loc[dt_quality_score['batch'] == 'QC_test', 'rank'].iloc[0])
        bad = "%.2f%s" % (len_bad/total_len * 100, '%')
        fair = "%.2f%s" % (len_fair/total_len * 100, '%')
        good = "%.2f%s" % (len_good/total_len * 100, '%')
        great = "%.2f%s" % (len_great/total_len * 100, '%')
        if query_rank == 1:
            queried = "%.2f%s" % (0, '%')
        elif query_rank == total_len + 1:
            queried = "%.2f%s" % (200, '%')
        else:
            queried = "%.2f%s" % (int(query_rank)/total_len *200, '%')
        
        # ticks number
        snr = dt_quality_score.loc[dt_quality_score['batch'] == 'QC_test', 'scaled_score'].iloc[0]
        Q0 = "%.1f" % float(dt_quality_score.loc[dt_quality_score['performance'] == 'Bad', 'scaled_score'].iloc[len_bad -1])
        Q1 = "%.1f" % float(dt_quality_score.loc[dt_quality_score['performance'] == 'Bad', 'scaled_score'].iloc[0])
        Q2 = "%.1f" % float(dt_quality_score.loc[dt_quality_score['performance'] == 'Fair', 'scaled_score'].iloc[0])
        Q3 = "%.1f" % float(dt_quality_score.loc[dt_quality_score['performance'] == 'Good', 'scaled_score'].iloc[0])
        Q4 = "%.1f" % float(dt_quality_score.loc[dt_quality_score['performance'] == 'Great', 'scaled_score'].iloc[0])
        
        # Position of ticks
        tick_Q1 = "%.2f%s" % (len_bad/total_len * 100, '%')
        tick_Q2 = "%.2f%s" % ((len_bad + len_fair)/total_len * 100, '%')
        tick_Q3 = "%.2f%s" % ((len_bad + len_fair +  len_good)/total_len * 100, '%')
        
        metrics_summary_html = """
        <!-- Arrow -->
        <div class="arrow" style="width: {queried}; margin-top:10px; height: 35px;">
        <svg class="lower-tangle" transform="translate(0, 18)"></svg>
        <span class="lower-label" style="margin-bottom: 25px;"><b> {snr} </b></span>
        </div>
        
        <!-- Progress bar -->
        <div class="progress">
          <div class="progress-bar progress-bar-bad" style="width: {bad}" data-toggle="tooltip" title="" data-original-title="">Bad</div>
          <div class="progress-bar progress-bar-fair" style="width: {fair}" data-toggle="tooltip" title="" data-original-title="">Fair</div>
          <div class="progress-bar progress-bar-good" style="width: {good}" data-toggle="tooltip" title="" data-original-title="">Good</div>
          <div class="progress-bar progress-bar-great" style="width: {great}" data-toggle="tooltip" title="" data-original-title="">Great</div>
        </div>
        
        <!-- Scale interval -->
        <span style="float:left; left:0%; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q0}</span>
        <span style="float:left; left:{tick_Q1}; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q1}</span>
        <span style="float:left; left:{tick_Q2}; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q2}</span>
        <span style="float:left; left:{tick_Q3}; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q3}</span>
        <span style="float:left; left:99%; position:relative; margin-top:-20px; color: #9F9FA3; font-size: 14px; text-align: center; display: inline-block">{Q4}</span>
        <br>
        """.format(bad=bad, fair=fair, good=good, great=great, queried=queried, snr=snr, tick_Q1=tick_Q1, tick_Q2=tick_Q2, tick_Q3=tick_Q3, Q0=Q0, Q1=Q1, Q2=Q2, Q3=Q3, Q4=Q4)

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
        'format': '{0:.3f}'
        }
        headers['historical_value'] = {
            'title': 'Historical value (mean ± SD)',
            'description': 'Historical Value (mean ± SD)',
            'scale': False,
            'format': '{0:.3f}'
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
            'col1_header': 'Quality Metrics',
            'no_beeswarm': True,
            'sortRows': False
            }
        
        metrics_table_html = table.plot(qc_metrics_summary_list, headers, table_config)

        # Add a report section with the scatter plot
        self.add_section(
            name='QC metrics summary',
            anchor=id + '_anchor',
            description="""
            The total performance score is calculated to measure the overall quality of a dataset generated from a lab for its effectiveness in quantifying the transcriptomic differences among the four Quartet RNA reference materials by summarizing reference dataset-independent quality measurement (SNR) and reference dataset-dependent quality measurement (RC). The total score is expressed as the geometrical mean of SNR and RC. 
            The submitted data to be tested can be divided into 4 levels based on the total score by comparing with historical batches: <span style="color: #b80d0d;font-weight:bold">Bad</span>, <span style="color: #d97c11;font-weight:bold">Fair</span>, <span style="color: #70c402;font-weight:bold">Good</span>, <span style="color: #0f9115;font-weight:bold">Great</span>.<br>
            * _Bad_ - the bottom 20%.
            * _Fair_ - between bottom 20% and bottom 50%.
            * _Good_ - between top 50% and top 20%.
            * _Great_ - the top 20%.
            """,
            plot = metrics_summary_html + '\n' +  metrics_table_html
            )

    # Plot2: snr and rc
    def plot_snr_rc_point(
            self,
            id,
            quality_score,
            title="SNR and RC",
            description=None,
            helptext=None):
        fig = px.scatter(quality_score, x="SNR", y="RC",
                         symbol = 'group', 
                         symbol_map = {"Reference": 0, "Query": 18},
                         color="group",
                         color_discrete_map={
                             "Query": "#bb1616",
                             "Reference": "#2f5c85"
                         },
                         marginal_x="box",
                         marginal_y="box",
                         template="simple_white")
        
        fig.update_traces(marker=dict(size=14))
        fig.update_layout(xaxis_title='SNR',
                          yaxis_title='RC',
                          font=dict(family="Arial, sans-serif",
                                    size=18,
                                    color="black"),
                          template="simple_white")
        fig.update_layout(legend_traceorder="reversed")
        fig.update_layout(legend_title_text='')

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='SNR and RC performance',
            anchor=id + '_anchor',
            description=description if description else
            'Performance metrics and thresholds using Quartet RNA reference materials',
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
        fig.update_layout(legend_title_text='')
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
            'Signal-to-noise ratio (SNR) is defined as the ratio of the power of a signal to the power of noise',
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
    
    # Plot4: RC correlation plot
    def plot_rc_cor_point(
            self,
            id,
            df_rc,
            title=None,
            description=None,
            helptext=None):
        cor_value = 'Correlation = ' + df_rc.iloc[1].at['cor'] + ' (N = ' + df_rc.iloc[1].at['gene_num'] + ')'
        fig = px.scatter(df_rc,
                         x="meanlogFC_ref",
                         y="meanlogFC_test",
                         color="compare",
                         color_discrete_map={
                             "D5/D6": "#4CC3D9",
                             "F7/D6": "#FFC65D",
                             "M8/D6": "#F16745"
                         },
                         marginal_x="box",
                         marginal_y="box",
                         template="simple_white")
        fig.update_layout(
            title=cor_value,
            xaxis_title='Reference Datasets',
            yaxis_title='Queried Data',
            font=dict(family="Arial, sans-serif", size=18, color="black"),
            template="simple_white")
        fig.update_layout(title={
            'y': 0.99,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        })
        fig.update_layout(legend_title_text='')

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='RC',
            anchor=id + '_anchor',
            description=description if description else
            'RC was calculate based on the Pearson correlation coefficient between the relative expression levels of a dataset for a given pair of groups and the corresponding reference fold-change values',
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
