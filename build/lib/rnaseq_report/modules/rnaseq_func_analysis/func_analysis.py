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
            name='功能预测',
            target="gene_expression_report",
            anchor='gene_expression_report',
            href='https://github.com/clinico-omics/rnaseq-report',
            info=" is an report module to show the plots and tables about functional analysis."
        )

        # Find and load any input files for DEG_GENE_PCA
        deg_kegg_data = []
        for f in self.find_log_files('rnaseq_deg_report/deg_kegg'):
            lines = f['f'].splitlines()
            keys = lines[0].split('\t')
            content = lines[1:]
            for values in content:
                deg_kegg_data.append(dict(zip(keys, values.split('\t'))))

        if len(deg_kegg_data) != 0:
            # Now add each section in order
            self.plot_deg_kegg('deg_kegg', deg_kegg_data)
        else:
            log.debug("No file matched: deg_kegg - kegg.txt")

        # Find and load any input files for DEG_GENE_VOLCANO
        deg_go_bp_data = []
        for f in self.find_log_files('rnaseq_deg_report/deg_go_bp'):
            deg_go_bp_data = f['f'].splitlines()

        if len(deg_go_bp_data) != 0:
            # Now add each section in order
            deg_go_bp_df = self.list2df(deg_go_bp_data)
            self.plot_deg_go_bp('deg_go_bp', deg_go_bp_df)
        else:
            log.debug("No file matched: deg_go_bp - go_bp.txt")

        # Find and load any input files for DEG_GENE_VOLCANO
        deg_go_cc_data = []
        for f in self.find_log_files('rnaseq_deg_report/deg_go_cc'):
            deg_go_cc_data = f['f'].splitlines()

        if len(deg_go_cc_data) != 0:
            # Now add each section in order
            deg_go_cc_df = self.list2df(deg_go_cc_data)
            self.plot_deg_go_cc('deg_go_cc', deg_go_cc_df)
        else:
            log.debug("No file matched: deg_go_cc - go_cc.txt")

        # Find and load any input files for DEG_GENE_VOLCANO
        deg_go_mf_data = []
        for f in self.find_log_files('rnaseq_deg_report/deg_go_mf'):
            deg_go_mf_data = f['f'].splitlines()

        if len(deg_go_mf_data) != 0:
            # Now add each section in order
            deg_go_mf_df = self.list2df(deg_go_mf_data)
            self.plot_deg_go_mf('deg_go_mf', deg_go_mf_df)
        else:
            log.debug("No file matched: deg_go_mf - go_mf.txt")

    def list2df(self, data):
        """Convert string list to dataframe"""
        cols = data[0].split('\t')
        array = []
        for line in data[1:]:
            array.append(line.split('\t'))
        df = pd.DataFrame(np.array(array), columns=cols)
        return df

    def plot_deg_kegg(self, id, deg_kegg_data, title="KEGG",
                      section_name="KEGG", description=None, helptext=None):
        """ Create the HTML for the deg pca plot """
        fig = px.bar(deg_kegg_data[0:19], x='Description', y='Count', color='p.adjust')

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })

        # Add a report section with the line plot
        self.add_section(
            name=section_name if section_name else 'KEGG',
            anchor=id + '_anchor',
            description=description if description else '分析差异表达基因在某一通路上是否过出现（over-presentation），横坐标为富集到的通路名称，纵坐标为富集到相应通路的基因数量，柱子的色阶代表相应通路富集显著性，颜色越深代表富集显著性越可靠。',
            helptext=helptext if helptext else '''
            富集是指将基因按照先验知识，也就是基因组注释信息，对基因进行分类的过程。基因经过分类后，能够帮助认知寻找到的基因是否具有某方面的共性(如功能、组成等等)。
            
            基因组的注释信息主要来自于数据库，如KEGG (yoto encyclopedia of genes and genomes, KEGG)和GO (genome annotation)等。

            KEGG主要的特点是将基因与各种生化反应联系在了一起。它提供的整合代谢途径查询十分出色，KEGG目前共包含了19个子数据库，他们被分类为系统信息、基因组信息和化学信息三个类别[6]。
            ''',
            plot=html
        )

    def plot_deg_go_bp(self, id, deg_go_bp_df, title="DEG GO BP",
                       section_name=None, description=None, helptext=None):
        fig = px.bar(deg_go_bp_df[0:19], x='Description', y='Count', color='p.adjust')

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })

        # Add a report section with the line plot
        self.add_section(
            name=section_name if section_name else 'GO BP',
            anchor=id + '_anchor',
            description=description if description else '分析差异表达基因在生物学过程（BP）中的富集分析，横坐标为富集到的细胞组分名称，纵坐标为富集到相应组分的基因数量，柱子的色阶代表相应通路富集显著性，颜色越深代表富集显著性越可靠。',
            helptext=helptext if helptext else '''
            基因本体数据库 (Gene Ontology) 是GO组织 (Gene Ontology Consortium) 构建的一个结构化的标准生物模型，旨在建立基因及其产物知识的标准词汇体系，涵盖了基因的细胞组分 (cellular component, CC)、分子功能 (molecular function, MF)、生物学过程(biological process, BP)[7]。
            
            Term是GO里面的基本描述单元，细胞组分（CC）指基因产物位于何种细胞器或基因产物组中(如糙面内质网，核糖体，蛋白酶体等)，即基因产物在什么地方起作用, 分子功能（MF）描述在个体分子生物学上的活性，如催化活性或结合活性，生物学过程（BP）由分子功能有序地组成的，具有多个步骤的一个过程。
            ''',
            plot=html
        )

    def plot_deg_go_cc(self, id, deg_go_cc_df, title="DEG GO CC",
                       section_name=None, description=None, helptext=None):
        fig = px.bar(deg_go_cc_df[0:19], x='Description', y='Count', color='p.adjust')

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })

        # Add a report section with the line plot
        self.add_section(
            name=section_name if section_name else 'GO CC',
            anchor=id + '_anchor',
            description=description if description else '分析差异表达基因在细胞组分（CC）中的富集分析，横坐标为富集到的细胞组分名称，纵坐标为富集到相应组分的基因数量，柱子的色阶代表相应通路富集显著性，颜色越深代表富集显著性越可靠。',
            helptext=helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot=html
        )

    def plot_deg_go_mf(self, id, deg_go_mf_df, title="DEG GO CC",
                       section_name=None, description=None, helptext=None):
        fig = px.bar(deg_go_mf_df[0:19], x='Description', y='Count', color='p.adjust')

        html = plotly_plot(fig, {
            'id': id + '_plot',
            'data_id': id + '_data',
            'title': title,
            'auto_margin': True
        })

        # Add a report section with the line plot
        self.add_section(
            name=section_name if section_name else 'GO MF',
            anchor=id + '_anchor',
            description=description if description else '分析差异表达基因在分子功能 （MF）中的富集分析，横坐标为富集到的细胞组分名称，纵坐标为富集到相应组分的基因数量，柱子的色阶代表相应通路富集显著性，颜色越深代表富集显著性越可靠。',
            helptext=helptext if helptext else '''
            This longer description explains what exactly the numbers mean
            and supports markdown formatting. This means that we can do _this_:

            * Something important
            * Something else important
            * Best of all - some `code`

            Doesn't matter if this is copied from documentation - makes it
            easier for people to find quickly.
            ''',
            plot=html
        )
