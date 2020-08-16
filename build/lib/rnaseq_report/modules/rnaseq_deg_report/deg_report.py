#!/usr/bin/env python

""" RnaSeqReport plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff
from .heatmap import heatmap

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
            name='差异表达分析',
            target="deg_report",
            anchor='deg_report',
            href='https://github.com/clinico-omics/rnaseq-report',
            info=" is an report module to show the plots and tables about differential expression gene."
        )

        # Find and load any input files for DEG_GENE_PCA
        deg_gene_pca_data = []
        for f in self.find_log_files('rnaseq_deg_report/deg_gene_pca'):
            lines = f['f'].splitlines()
            keys = lines[0].split('\t')
            content = lines[1:]
            for values in content:
                deg_gene_pca_data.append(dict(zip(keys, values.split('\t'))))

        if len(deg_gene_pca_data) != 0:
            ## Now add each section in order
            self.plot_deg_gene_pca('deg_pca_plot', deg_gene_pca_data)
        else:
            log.debug("No file matched: deg_gene_pca - deg_pca.txt")

        # Find and load any input files for DEG_GENE_VOLCANO
        deg_gene_volcano_data = []
        for f in self.find_log_files('rnaseq_deg_report/deg_gene_volcano'):
            deg_gene_volcano_data = f['f'].splitlines()

        if len(deg_gene_volcano_data) != 0:
            ## Now add each section in order
            deg_gene_volcano_df = self.list2df(deg_gene_volcano_data)
            self.plot_deg_gene_volcano('deg_gene_volcano', deg_gene_volcano_df)
        else:
            log.debug("No file matched: deg_gene_volcano - deg_fc.txt")

        # Find and load any input files for DEG_GENE_HEATMAP
        deg_gene_heatmap_data = []
        for f in self.find_log_files('rnaseq_deg_report/deg_gene_heatmap'):
            deg_gene_heatmap_data = f['f'].splitlines()

        if len(deg_gene_heatmap_data) != 0:
            ## Now add each section in order
            deg_gene_heatmap_df = self.list2df(deg_gene_heatmap_data, drop=[0,])
            self.plot_deg_gene_heatmap('deg_gene_heatmap', deg_gene_heatmap_df)
        else:
            log.debug("No file matched: deg_gene_heatmap - deg_exp.txt")

    def remove_by_idx(self, data, drop=[]):
        for idx in drop:
            del data[idx]
        return data

    def list2df(self, data, drop=[]):
        """Convert string list to dataframe"""
        cols = self.remove_by_idx(data[0].split('\t'), drop)
        array = []
        for line in data[1:]:
            if len(drop) != 0:
                array.append(self.remove_by_idx(line.split('\t'), drop))
            else:
                array.append(line.split('\t'))
        df = pd.DataFrame(np.array(array), columns=cols)
        return df

    def plot_deg_gene_pca(self, id, deg_gene_pca_data, title="DEG PCA Plot",
                          section_name="DEG PCA Plot", description=None, helptext=None):
        """ Create the HTML for the deg pca plot """
        color_palette = ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f']

        groups = list(map(lambda x: x['group'], deg_gene_pca_data))
        if len(groups) > 6:
            log.error('Groups is greater than 6.')
            colors = None
        else:
            colors = dict(zip(groups, color_palette))

        data_labels = list(map(lambda x: {'name': x['group']}, deg_gene_pca_data))

        pconfig = {
            'id': id + '_plot',
            'title': title,
            'ylab': 'PCA2',
            'xlab': 'PCA1',
            'tt_label': 'PC1 {point.x:.2f}<br>PC2 {point.y:.2f}',
            'data_labels': data_labels
        }

        data = dict()
        for item in deg_gene_pca_data:
            group_name = item['group']
            point = {
                'x': float(item['PC1']),
                'y': float(item['PC2']),
                'name': item['SAMPLE_ID'],
                'color': colors.get(group_name) if colors else '#a6cee3'
            }

            if data.get(group_name):
                data[group_name].append(point)
            else:
                data[group_name] = [point, ]

        scatter_plot_html = scatter.plot(data, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name = 'PCA',
            anchor = id + '_anchor',
            description = description if description else '差异表达基因的表达矩阵进行主成分分析（PCA），通过降维的方式将基因表达矩阵中在更小的维度下展示数据的特征，把一系列可能线性相关的变量转换为一组线性不相关的新变量，也称为主成分，利用主成分分析样品之间相关性，确定样品总体上的差异，或者查看是否有批次效应。',
            helptext = helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot = scatter_plot_html
        )

    def plot_deg_gene_volcano(self, id, deg_gene_volcano_df, title="DEG Gene Volcano",
                              section_name=None, description=None, helptext=None):
        fig = px.scatter(deg_gene_volcano_df, x="logFC", y="PvalueLog", color="group")

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })

        # Add a report section with the line plot
        self.add_section(
            name = section_name if section_name else '火山图',
            anchor = id + '_anchor',
            description = description if description else '火山图可以方便直观地展示两组样本间基因差异表达的分布情况，通常横坐标用log2（fold change）表示，差异越大的基因分布在两端，纵坐标用-log10（pvalue）表示，T检验显著性P值的负对数，通常差异倍数越大的基因T检验越显著，所以往往关注左上角和右上角的值。',
            helptext = helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot = html
        )

    def plot_deg_gene_heatmap(self, id, deg_gene_heatmap_df, title="DEG Gene Heatmap",
                              section_name=None, description=None, helptext=None):
        fig = heatmap(deg_gene_heatmap_df)

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })

        # Add a report section with the line plot
        self.add_section(
            name = section_name if section_name else 'DEG Gene Heatmap',
            anchor = id + '_anchor',
            description = description if description else 'This plot shows some numbers, and how they relate.',
            helptext = helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot = html
        )