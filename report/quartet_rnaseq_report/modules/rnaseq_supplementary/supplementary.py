""" Quartet DNAseq Report plugin module """

from __future__ import print_function
import logging
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='Supplementary',
            target='supplementary',
            anchor='supplementary',
            href='https://github.com/clinico-omics/quartet-rnaseq-report',
            info=
            ' is a module to show the additional information about this quality assessment report.'
        )

        html = '''
            <!-- Method -->
            <div class='method'>
                <div class='small-12 columns'>
                    <h3 class='section-header black'>Method</h3>
                    <img src='quartet_dnaseq_report/modules/supplementary/assets/img/overall_process_of_data_analysis' title='overall process of data analysis' width='70%' height='70%'/>
                    <p>
                        We use FastQC, FastQ Screen, Qualimap and MultiQC to evaluate the quality of sequencing data. [<a class='reference' href='#ref-1'>1</a>] RNA-seq quality control consists of pre-alignment, post-alignment and quantification quality control.
                    </p>
                    <p>
                        Pre-alignment quality control focuses on raw fastq files and helps to determine systematic bias and library issue, such as sequencing quality issue, high GC or AT, PCR bias, adapter contaminant, cross species contamination. Fastqc [<a class='reference' href='#ref-2'>2</a>] and fastqscreen [<a class='reference' href='#ref-3'>3</a>] are used to evaluate raw reads quality.
                    </p>
                    <p>
                        Post-alignment quality control focuses on bam files and helps to measure library performance and sample variance, such as sequencing error rate, sequencing depth and coverage consistency. Qualimap [<a class='reference' href='#ref-4'>4</a>] is used to evaluate quality of bam files.
                    </p>
                    <p>
                        Quantification quality control is to evaluate the data quality from a set of data quality control metrics and thresholds, especially discrimination of different groups. 
                    </p>
                </div>
            </div>
            <!-- Reference -->
            <div class='reference'>
                <div class='small-12 columns'>
                <h3 class='section-header black'>References</h3>
                <p class='reference'><a name='ref-1'>1.</a> <a href='https://www.bioinformatics.babraham.ac.uk/projects/fastqc/'>https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a></p>
                <p class='reference'><a name='ref-2'>2.</a> <a href='https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/'>https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/</a></p>
                <p class='reference'><a name='ref-3'>3.</a> <a href='http://qualimap.bioinfo.cipf.es/'>http://qualimap.bioinfo.cipf.es/</a></p>
                <p class='reference'><a name='ref-4'></a>4. Ewels P, et al. Bioinformatics, 2016. </p>
                <p class='reference'><a name='ref-5'></a>5. Pertea M, et.al. Nature Biotechnology, 2015 </p>
                </div>
            </div>
            <!-- Software version -->
            <div class='software'>
                <h3 class='section-header black'>Software</h3>
                <dl class='dl-horizontal'>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>Fastqc</dt><dd>v0.11.5</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>Fastqscreen</dt><dd>v0.12.0</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>Qualimap</dt><dd>v2.0.0</dd>
                    <dt style='text-align:left;width:300px;font-weight:normal;margin-top:1ex'>MultiQC</dt><dd>v1.9</dd>
                </dl>
            </div>
            <!-- Contact us -->
            <div class='contact'>
                <div class='small-12 columns'>
                <h3 class='section-header black'>Contact Us</h3>
                <b>Fudan University Pharmacogenomics Research Center</b>
                <p><strong>Project Manager Zhihui Li</strong></p>
                <li style='margin-top:1ex'>Phone: 15200852771</li>
                <li style='margin-top:1ex'>Email: 18210700119@fudan.edu.cn</li>
                </div>
            </div>
            <!-- Disclaimer -->
            <div class='disclaimer'>
                <div class='small-12 columns'>
                <h3 class='section-header black'>Disclaimer</h3>
                <p>This quality control report is only for this specific test data set and doesn’t represent an evaluation of the business level of the sequencing company. This report is only used for scientific research, not for clinical or commercial use. We don’t bear any economic and legal liabilities for any benefits or losses (direct or indirect) from using the results of this report.</p>
                </div>
            </div>
            '''

        self.add_section(name='', anchor='', description='', plot=html)