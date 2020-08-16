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
            name='Performance Assessment',
            target="performance_assessment",
            anchor='performance_assessment',
            href='https://github.com/clinico-omics/rnaseq-report',
            info=
            " is an report module to show the performance of quatet samples.")

        # Find and load any input files for DEG_GENE_PCA
        reference_deg_plot_data = []
        for f in self.find_log_files(
                'performance_assessment/ref_degs_performance_compared'):
            lines = f['f'].splitlines()
            keys = lines[0].split('\t')
            content = lines[1:]
            for values in content:
                ref_degs_data.append(dict(zip(keys, values.split('\t'))))

        if len(ref_degs_data) != 0:
            ## Now add each section in order
            self.plot_ref_degs_point('ref_degs_plot', ref_degs_data)
        else:
            log.debug(
                "No file matched: ref_degs_point - ref_degs_performance_compared.txt"
            )

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

    def plot_ref_degs_point(self,
                            id,
                            reference_deg_plot_data,
                            title="Reference DEGs Performance Compared",
                            section_name=None,
                            description=None,
                            helptext=None):
        fig = px.scatter(ref_degs_data,
                         x="Sensitivity",
                         y="Specificity",
                         color="Group")

        html = plotly_plot(
            fig, {
                'id': id + '_plot',
                'data_id': id + '_data',
                'title': title,
                'auto_margin': True
            })

        # Add a report section with the line plot
        self.add_section(
            name='Performance Assessment',
            anchor=id + '_anchor',
            description=description if description else
            '差异表达基因的表达矩阵进行主成分分析（PCA），通过降维的方式将基因表达矩阵中在更小的维度下展示数据的特征，把一系列可能线性相关的变量转换为一组线性不相关的新变量，也称为主成分，利用主成分分析样品之间相关性，确定样品总体上的差异，或者查看是否有批次效应。',
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
