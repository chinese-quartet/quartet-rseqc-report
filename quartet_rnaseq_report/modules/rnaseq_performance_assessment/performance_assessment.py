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
from rnaseq_report.modules.plotly import plot as plotly_plot

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Performance Assessment',
            target="performance_assessment",
            anchor='performance_assessment',
            href='https://github.com/clinico-omics/rnaseq-report',
            info=
            " is an report module to show the performance of quartet samples.")
        ########################find log files#################################
        # Plot0: Performance Score
        performance_score = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/performance_score'):
            performance_score = f['f'].splitlines()

        if len(performance_score) != 0:
            ## Now add each section in order
            performance_score_df = self.list2df(performance_score)
            test_score = performance_score_df.loc[
                performance_score_df['Batch'] ==
                'Test']['Benchmarkingscore'].values[0]
            performance_score_list = [
                performance_score_df['Benchmarkingscore'].values.tolist()
            ]
            self.plot_performance_score('plot_performance_score',
                                        performance_score_list, test_score)
        else:
            log.debug(
                "No file matched: performance_score - performance_score.txt")

        # Summary Table1: DEG performance summary table
        deg_performance_summary = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/deg_performance_summary'):
            deg_performance_summary = f['f'].splitlines()

        if len(deg_performance_summary) != 0:
            ## Now add each section in order
            deg_performance_summary_df = self.list2df(deg_performance_summary)
            deg_performance_summary_list = [
                deg_performance_summary_df.columns.values.tolist()
            ] + deg_performance_summary_df.values.tolist()
            self.plot_deg_summary_table('deg_summary_table',
                                        deg_performance_summary_list)
        else:
            log.debug(
                "No file matched: deg_performance_summary_table - deg_performance_summary.txt"
            )

        # Plot1: Relative expression performance evaluation
        ref_degs_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/ref_degs_performance_compared'):
            ref_degs_data = f['f'].splitlines()

        if len(ref_degs_data) != 0:
            ## Now add each section in order
            ref_degs_data_df = self.list2df(ref_degs_data)
            self.plot_ref_degs_point('ref_degs_plot', ref_degs_data_df,
                                     deg_performance_summary_list)
        else:
            log.debug(
                "No file matched: ref_degs_point - ref_degs_performance_compared.txt"
            )

        # Summary Table2: Relative expression performance summary table
        rel_performance_summary = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/rel_performance_summary'):
            rel_performance_summary = f['f'].splitlines()

        if len(rel_performance_summary) != 0:
            ## Now add each section in order
            rel_performance_summary_df = self.list2df(rel_performance_summary)
            rel_performance_summary_list = [
                rel_performance_summary_df.columns.values.tolist()
            ] + rel_performance_summary_df.values.tolist()
            self.plot_rel_summary_table('rel_summary_table',
                                        rel_performance_summary_list)
        else:
            log.debug(
                "No file matched: rel_performance_summary_table - rel_performance_summary.txt"
            )

        # Plot2: Relative expression performance evaluation
        ref_rel_exp_per_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/ref_rel_exp_per_compared'):
            ref_rel_exp_per_data = f['f'].splitlines()

        if len(ref_rel_exp_per_data) != 0:
            ## Now add each section in order
            ref_rel_exp_per_data_df = self.list2df(ref_rel_exp_per_data)
            self.plot_ref_rel_exp_per_point('ref_rel_exp_per_plot',
                                            ref_rel_exp_per_data_df)
        else:
            log.debug(
                "No file matched: plot_ref_rel_exp_per_point - ref_rel_exp_per_compared.txt"
            )

        # Summary Table3: Study design performance summary table
        studydesign_performance_summary = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/studydesign_performance_summary'
        ):
            studydesign_performance_summary = f['f'].splitlines()

        if len(studydesign_performance_summary) != 0:
            ## Now add each section in order
            studydesign_performance_summary_df = self.list2df(
                studydesign_performance_summary)
            studydesign_performance_summary_list = [
                studydesign_performance_summary_df.columns.values.tolist()
            ] + studydesign_performance_summary_df.values.tolist()
            self.plot_studydesign_summary_table(
                'studydesign_summary_table',
                studydesign_performance_summary_list)
        else:
            log.debug(
                "No file matched: studydesign_performance_summary_table - studydesign_performance_summary.txt"
            )

        # Plot3: SNR
        snr_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/studydesign_snr'):
            snr_data = f['f'].splitlines()

        if len(snr_data) != 0:
            ## Now add each section in order
            snr_data_df = self.list2df(snr_data)
            self.plot_snr_pca_point('snr_pca_plot', snr_data_df)
        else:
            log.debug("No file matched: snr_pca_point - studydesign_snr.txt")

        # Plot4: D5 correlation
        d5_cor_data = []
        for f in self.find_log_files(
                'rnaseq_performance_assessment/d5_correlation'):
            d5_cor_data = f['f'].splitlines()

        if len(d5_cor_data) != 0:
            ## Now add each section in order
            d5_cor_data_df = self.list2df(d5_cor_data)
            self.plot_d5_cor_point('d5_cor_plot', d5_cor_data_df)
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
    ###Plot0: performance_score
    def plot_performance_score(self,
                               id,
                               performance_score_list,
                               test_score,
                               title=None,
                               section_name=None,
                               description=None,
                               helptext=None):
        fig = px.imshow(performance_score_list,
                        x=performance_score_list[0],
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
            name='Performance Score',
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

    ###Table1: DEG performance summary
    def plot_deg_summary_table(self,
                               id,
                               deg_performance_summary_list,
                               title="DEG Performance Summary",
                               section_name=None,
                               description=None,
                               helptext=None):
        fig = ff.create_table(deg_performance_summary_list, height_constant=60)

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='DEG performance summary',
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

    # Plot1: DEGs performance evaluation
    def plot_ref_degs_point(self,
                            id,
                            ref_degs_data_df,
                            deg_performance_summary_list,
                            title="Reference DEG Performance Compared",
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
            name='DEG Performance Based On Reference',
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

    # Table2: Relative expression performance summary
    def plot_rel_summary_table(self,
                               id,
                               rel_performance_summary_list,
                               title="Relative Expression Performance Summary",
                               section_name=None,
                               description=None,
                               helptext=None):
        fig = ff.create_table(rel_performance_summary_list, height_constant=60)

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='Relative Expression Performance Summary',
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
    def plot_ref_rel_exp_per_point(self,
                                   id,
                                   ref_rel_exp_per_data_df,
                                   title="Reference DEGs Performance Compared",
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

    # Table3: Study design performance summary
    def plot_studydesign_summary_table(
            self,
            id,
            studydesign_performance_summary_list,
            title="Study Design Performance Summary",
            section_name=None,
            description=None,
            helptext=None):
        fig = ff.create_table(studydesign_performance_summary_list,
                              height_constant=60)

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the scatter plot
        self.add_section(
            name='Study Design Performance Summary',
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
                           title="SNR",
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
            name='SNR',
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

    # Plot4: Correlation figure
    def plot_d5_cor_point(self,
                          id,
                          d5_cor_data_df,
                          title="Correlation of expression",
                          section_name=None,
                          description=None,
                          helptext=None):
        fig = px.scatter(d5_cor_data_df,
                         x="D5_1",
                         y="D5_2",
                         hover_name='GENE_ID',
                         template="simple_white")

        fig.update_layout(title='CTR',
                          xaxis_title='D5_1',
                          yaxis_title='D5_2',
                          font=dict(family="Arial, sans-serif",
                                    size=18,
                                    color="black"))

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
            name='Correlation of expression',
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
