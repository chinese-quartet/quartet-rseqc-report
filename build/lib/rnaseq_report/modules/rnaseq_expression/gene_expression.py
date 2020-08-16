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
            name='基因表达',
            target="gene_expression_report",
            anchor='gene_expression_report',
            href='https://github.com/clinico-omics/rnaseq-report',
            info=" 基因表达即基因从DNA到RNA转录和蛋白翻译的过程，转录本的丰度体现基因的表达水平，转录本丰度越高，则基因表达水平越高。在RNA-seq分析中，通过计算比对到基因组特定区域的序列(clean reads)来衡量基因或转录本的表达水平。"
        )

        # Find and load any input files for DEG_GENE_PCA
        all_gene_boxplot_data = []
        for f in self.find_log_files('rnaseq_deg_report/all_gene_boxplot'):
            lines = f['f'].splitlines()
            keys = lines[0].split('\t')
            content = lines[1:]
            for values in content:
                all_gene_boxplot_data.append(dict(zip(keys, values.split('\t'))))

        if len(all_gene_boxplot_data) != 0:
            ## Now add each section in order
            self.plot_all_gene_boxplot('all_gene_boxplot', all_gene_boxplot_data)
        else:
            log.debug("No file matched: all_gene_boxplot - all_gene_exp.txt")

        # Find and load any input files for DEG_GENE_VOLCANO
        all_gene_pca_data = []
        for f in self.find_log_files('rnaseq_deg_report/all_gene_pca'):
            all_gene_pca_data = f['f'].splitlines()

        if len(all_gene_pca_data) != 0:
            ## Now add each section in order
            all_gene_pca_df = self.list2df(all_gene_pca_data)
            self.plot_all_gene_pca('all_gene_pca', all_gene_pca_df)
        else:
            log.debug("No file matched: all_gene_pca - all_gene_pca.txt")

    def list2df(self, data):
        """Convert string list to dataframe"""
        cols = data[0].split('\t')
        array = []
        for line in data[1:]:
            array.append(line.split('\t'))
        df = pd.DataFrame(np.array(array), columns=cols)
        return df

    def plot_all_gene_boxplot(self, id, all_gene_boxplot_data, title="基因表达分布",
                              section_name="基因表达分布", description=None, helptext=None):
        """ Create the HTML for the deg pca plot """
        fig = px.box(all_gene_boxplot_data, x="sample_id", y="expression", color="group", notched=True)

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })

        # Add a report section with the line plot
        self.add_section(
            name = section_name if section_name else '基因表达分布',
            anchor = id + '_anchor',
            description = description if description else '通过箱型图展示每个样本中基因表达分布，横轴是样本，纵轴是每个基因的表达，通过展示每个样本中基因表达分布来说明每个样本基因定量水平是否一致。',
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

    def plot_all_gene_pca(self, id, all_gene_pca_df, title="PCA",
                          section_name=None, description=None, helptext=None):
        fig = px.scatter(all_gene_pca_df, x="PC1", y="PC2", color="group")

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })

        # Add a report section with the line plot
        self.add_section(
            name = section_name if section_name else 'PCA',
            anchor = id + '_anchor',
            description = description if description else '所有基因的表达矩阵进行主成分分析（PCA），通过降维的方式将基因表达矩阵中在更小的维度下展示数据的特征，把一系列可能线性相关的变量转换为一组线性不相关的新变量，也称为主成分，利用主成分分析样品之间相关性，确定样品总体上的差异，或者查看是否有批次效应。',
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
