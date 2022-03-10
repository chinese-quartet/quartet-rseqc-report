#!/usr/bin/env python
""" RnaSeqReport plugin module """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import logging
import os
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ Qualimap is really a collection of separate programs:
    BamQC, RNASeq and Counts.. This module is split into separate
    files to reflect this and help with code organisation. """
    def __init__(self):

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Post Alignment QC',
            anchor='rnaseq_post_alignment_qc',
            href='https://github.com/clinico-omics/rnaseq-report',
            info=" is an report module to show the quality of base.")

        # Initialise the submodules
        from . import QM_BamQC

        # Set up class objects to hold parsed data()
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = defaultdict(lambda: dict())
        n = dict()

        n['BamQC'] = QM_BamQC.parse_reports(self)
        if n['BamQC'] > 0:
            log.info("Found {} BamQC reports".format(n['BamQC']))

        ### add report section-------------
        # Make the plots for the report
        self.add_qm_rnaseq_stats_data()
        self.post_alignment_stats()
        QM_BamQC.report_sections(self)
        self.qm_rnaseq_parse_reports()

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

    # Helper functions
    def get_s_name(self, f):
        s_name = os.path.basename(os.path.dirname(f['root']))
        s_name = self.clean_s_name(s_name, f['root'])
        if s_name.endswith('.qc'):
            s_name = s_name[:-3]
        return s_name

    def post_alignment_stats(self):
        """ Add some single-number stats to the basic statistics
        table at the top of the report """
        # Prep the data
        error_rate_mean = []
        gc_mean = []
        median_insert_size_mean = []
        exonic_mean = []
        intronic_mean = []
        intergenic_mean = []

        for k in self.general_stats_data.keys():
            error_rate_mean.append(
                self.general_stats_data[k]['general_error_rate'])
            gc_mean.append(self.general_stats_data[k]['avg_gc'])
            median_insert_size_mean.append(
                self.general_stats_data[k]['median_insert_size'])
            exonic_mean.append(
                self.general_stats_data[k]['reads_aligned_exonic'])
            intronic_mean.append(
                self.general_stats_data[k]['reads_aligned_intronic'])
            intergenic_mean.append(
                self.general_stats_data[k]['reads_aligned_intergenic'])

        self.general_stats_data['Batch average value'] = {
            'general_error_rate':
            sum(error_rate_mean) / len(error_rate_mean),
            'avg_gc':
            sum(gc_mean) / len(gc_mean),
            'median_insert_size':
            sum(median_insert_size_mean) / len(median_insert_size_mean),
            'reads_aligned_exonic':
            sum(exonic_mean) / len(exonic_mean),
            'reads_aligned_intronic':
            sum(intronic_mean) / len(intronic_mean),
            'reads_aligned_intergenic':
            sum(intergenic_mean) / len(intergenic_mean)
        }

        self.general_stats_data['Historical value'] = {
            'general_error_rate': '0.01±0.01'
        }

        self.add_section(name='Post alignment stats',
                         anchor='post_alignment_stats',
                         description='The summary of qualimap',
                         helptext='''
            To enable multiple samples to be plotted on the same graph, only the mean quality
            scores are plotted (unlike the box plots seen in FastQC reports).
            Taken from the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html):
            _The y-axis on the graph shows the quality scores. The higher the score, the better
            the base call. The background of the graph divides the y axis into very good quality
            calls (green), calls of reasonable quality (orange), and calls of poor quality (red).
            The quality of calls on most platforms will degrade as the run progresses, so it is
            common to see base calls falling into the orange area towards the end of a read._
            ''',
                         plot=table.plot(self.general_stats_data,
                                         self.general_stats_headers))

    def add_qm_rnaseq_stats_data(self):

        regexes = {
            'reads_aligned': r"read(?:s| pairs) aligned\s*=\s*([\d,]+)",
            'total_alignments': r"total alignments\s*=\s*([\d,]+)",
            'non_unique_alignments': r"non-unique alignments\s*=\s*([\d,]+)",
            'reads_aligned_genes': r"aligned to genes\s*=\s*([\d,]+)",
            'ambiguous_alignments': r"ambiguous alignments\s*=\s*([\d,]+)",
            'not_aligned': r"not aligned\s*=\s*([\d,]+)",
            '5_3_bias': r"5'-3' bias\s*=\s*([\d,\.]+)$",
            'reads_aligned_exonic': r"exonic\s*=\s*([\d,]+)",
            'reads_aligned_intronic': r"intronic\s*=\s*([\d,]+)",
            'reads_aligned_intergenic': r"intergenic\s*=\s*([\d,]+)",
            'reads_aligned_overlapping_exon':
            r"overlapping exon\s*=\s*([\d,]+)",
        }
        for f in self.find_log_files(
                'rnaseq_post_alignment_qc/rnaseq_qc/rnaseq_qc_results'):
            d = dict()

            # Get the sample name
            s_name_regex = re.search(r"bam file\s*=\s*(.+)", f['f'],
                                     re.MULTILINE)
            if s_name_regex:
                d['bam_file'] = s_name_regex.group(1)
                s_name = self.clean_s_name(d['bam_file'], f['root'])
            else:
                log.warn(
                    "Couldn't find an input filename in genome_results file {}/{}"
                    .format(f['root'], f['fn']))
                return None

            # Check for and 'fix' European style decimal places / thousand separators
            comma_regex = re.search(r"exonic\s*=\s*[\d\.]+ \(\d{1,3},\d+%\)",
                                    f['f'], re.MULTILINE)
            if comma_regex:
                log.debug(
                    "Trying to fix European comma style syntax in Qualimap report {}/{}"
                    .format(f['root'], f['fn']))
                f['f'] = f['f'].replace('.', '')
                f['f'] = f['f'].replace(',', '.')

            # Go through all numeric regexes
            for k, r in regexes.items():
                r_search = re.search(r, f['f'], re.MULTILINE)
                if r_search:
                    try:
                        d[k] = float(r_search.group(1).replace(',', ''))
                    except UnicodeEncodeError:
                        # Qualimap reports infinity (\u221e) when 3' bias denominator is zero
                        pass
                    except ValueError:
                        d[k] = r_search.group(1)

            for k in [
                    'reads_aligned_exonic', 'reads_aligned_intronic',
                    'reads_aligned_intergenic'
            ]:
                try:
                    self.general_stats_data[s_name][k] = d[k]
                except KeyError:
                    pass

        #### General Stats
        self.general_stats_headers['reads_aligned_exonic'] = {
            'title': 'Exonic',
            'description': 'Reads Aligned Exonic',
            'min': 0,
            'scale': 'RdBu',
            'shared_key': 'read_count'
        }

        self.general_stats_headers['reads_aligned_intronic'] = {
            'title': 'Intronic',
            'description': 'Reads Aligned Intronic',
            'min': 0,
            'scale': 'RdBu',
            'shared_key': 'read_count'
        }

        self.general_stats_headers['reads_aligned_intergenic'] = {
            'title': 'Intergenic',
            'description': 'Reads Aligned Intergenic',
            'min': 0,
            'scale': 'RdBu',
            'shared_key': 'read_count'
        }

    def qm_rnaseq_parse_reports(self):
        """ Find Qualimap RNASeq reports and parse their data """

        self.qualimap_rnaseq_genome_results = dict()
        regexes = {
            'reads_aligned': r"read(?:s| pairs) aligned\s*=\s*([\d,]+)",
            'total_alignments': r"total alignments\s*=\s*([\d,]+)",
            'non_unique_alignments': r"non-unique alignments\s*=\s*([\d,]+)",
            'reads_aligned_genes': r"aligned to genes\s*=\s*([\d,]+)",
            'ambiguous_alignments': r"ambiguous alignments\s*=\s*([\d,]+)",
            'not_aligned': r"not aligned\s*=\s*([\d,]+)",
            '5_3_bias': r"5'-3' bias\s*=\s*([\d,\.]+)$",
            'reads_aligned_exonic': r"exonic\s*=\s*([\d,]+)",
            'reads_aligned_intronic': r"intronic\s*=\s*([\d,]+)",
            'reads_aligned_intergenic': r"intergenic\s*=\s*([\d,]+)",
            'reads_aligned_overlapping_exon':
            r"overlapping exon\s*=\s*([\d,]+)",
        }
        for f in self.find_log_files(
                'rnaseq_post_alignment_qc/rnaseq_qc/rnaseq_qc_results'):
            d = dict()

            # Get the sample name
            s_name_regex = re.search(r"bam file\s*=\s*(.+)", f['f'],
                                     re.MULTILINE)
            if s_name_regex:
                d['bam_file'] = s_name_regex.group(1)
                s_name = self.clean_s_name(d['bam_file'], f['root'])
            else:
                log.warn(
                    "Couldn't find an input filename in genome_results file {}/{}"
                    .format(f['root'], f['fn']))
                return None

            # Check for and 'fix' European style decimal places / thousand separators
            comma_regex = re.search(r"exonic\s*=\s*[\d\.]+ \(\d{1,3},\d+%\)",
                                    f['f'], re.MULTILINE)
            if comma_regex:
                log.debug(
                    "Trying to fix European comma style syntax in Qualimap report {}/{}"
                    .format(f['root'], f['fn']))
                f['f'] = f['f'].replace('.', '')
                f['f'] = f['f'].replace(',', '.')

            # Go through all numeric regexes
            for k, r in regexes.items():
                r_search = re.search(r, f['f'], re.MULTILINE)
                if r_search:
                    try:
                        d[k] = float(r_search.group(1).replace(',', ''))
                    except UnicodeEncodeError:
                        # Qualimap reports infinity (\u221e) when 3' bias denominator is zero
                        pass
                    except ValueError:
                        d[k] = r_search.group(1)

            # Save results
            if s_name in self.qualimap_rnaseq_genome_results:
                log.debug(
                    "Duplicate genome results sample name found! Overwriting: {}"
                    .format(s_name))
            self.qualimap_rnaseq_genome_results[s_name] = d
            self.add_data_source(f,
                                 s_name=s_name,
                                 section='rna_genome_results')

        #### Coverage profile
        self.qualimap_rnaseq_cov_hist = dict()
        for f in self.find_log_files(
                'rnaseq_post_alignment_qc/rnaseq_qc/coverage',
                filehandles=True):
            s_name = self.get_s_name(f)
            d = dict()
            for l in f['f']:
                if l.startswith('#'):
                    continue
                coverage, count = l.split(None, 1)
                coverage = int(round(float(coverage)))
                count = float(count)
                d[coverage] = count

            if len(d) == 0:
                log.debug(
                    "Couldn't parse contents of coverage histogram file {}".
                    format(f['fn']))
                return None

            # Save results
            if s_name in self.qualimap_rnaseq_cov_hist:
                log.debug(
                    "Duplicate coverage histogram sample name found! Overwriting: {}"
                    .format(s_name))
            self.qualimap_rnaseq_cov_hist[s_name] = d
            self.add_data_source(f,
                                 s_name=s_name,
                                 section='rna_coverage_histogram')

        # Filter to strip out ignored sample names
        self.qualimap_rnaseq_genome_results = self.ignore_samples(
            self.qualimap_rnaseq_genome_results)
        self.qualimap_rnaseq_cov_hist = self.ignore_samples(
            self.qualimap_rnaseq_cov_hist)

        #### Plots

        # Genomic Origin Bar Graph
        # NB: Ignore 'Overlapping Exon' in report - these make the numbers add up to > 100%
        if len(self.qualimap_rnaseq_genome_results) > 0:
            gorigin_cats = OrderedDict()
            gorigin_cats['reads_aligned_exonic'] = {'name': 'Exonic'}
            gorigin_cats['reads_aligned_intronic'] = {'name': 'Intronic'}
            gorigin_cats['reads_aligned_intergenic'] = {'name': 'Intergenic'}
            gorigin_pconfig = {
                'id': 'qualimap_genomic_origin',
                'title': 'Qualimap RNAseq: Genomic Origin',
                'ylab': 'Number of reads',
                'cpswitch_c_active': False
            }
            genomic_origin_helptext = '''
            There are currently three main approaches to map reads to transcripts in an
            RNA-seq experiment: mapping reads to a reference genome to identify expressed
            transcripts that are annotated (and discover those that are unknown), mapping
            reads to a reference transcriptome, and <i>de novo</i> assembly of transcript
            sequences (<a href="https://doi.org/10.1186/s13059-016-0881-8"
            target="_blank">Conesa et al. 2016</a>).

            For RNA-seq QC analysis, QualiMap can be used to assess alignments produced by
            the first of these approaches. For input, it requires a GTF annotation file
            along with a reference genome, which can be used to reconstruct the exon
            structure of known transcripts. This allows mapped reads to be grouped by
            whether they originate in an exonic region (for QualiMap, this may include
            5&#8242; and 3&#8242; UTR regions as well as protein-coding exons), an intron,
            or an intergenic region (see the <a href="http://qualimap.bioinfo.cipf.es/doc_html/index.html"
            target="_blank">Qualimap 2 documentation</a>).

            The inferred genomic origins of RNA-seq reads are presented here as a bar graph
            showing either the number or percentage of mapped reads in each read dataset
            that have been assigned to each type of genomic region. This graph can be used
            to assess the proportion of useful reads in an RNA-seq experiment. That
            proportion can be reduced by the presence of intron sequences, especially if
            depletion of ribosomal RNA was used during sample preparation (<a href="https://doi.org/10.1038/nrg3642"
            target="_blank">Sims et al. 2014</a>). It can also be reduced by off-target
            transcripts, which are detected in greater numbers at the sequencing depths
            needed to detect poorly-expressed transcripts (<a href="https://doi.org/10.1101/gr.124321.111"
            target="_blank">Tarazona et al. 2011</a>).'''
            self.add_section(
                name='Genomic origin of reads',
                anchor='qualimap-reads-genomic-origin',
                description=
                'Classification of mapped reads as originating in exonic, intronic or intergenic regions. These can be displayed as either the number or percentage of mapped reads.',
                helptext=genomic_origin_helptext,
                plot=bargraph.plot(self.qualimap_rnaseq_genome_results,
                                   gorigin_cats, gorigin_pconfig))

        if len(self.qualimap_rnaseq_cov_hist) > 0:
            coverage_profile_helptext = '''
            There are currently three main approaches to map reads to transcripts in an
            RNA-seq experiment: mapping reads to a reference genome to identify expressed
            transcripts that are annotated (and discover those that are unknown), mapping
            reads to a reference transcriptome, and <i>de novo</i> assembly of transcript
            sequences (<a href="https://doi.org/10.1186/s13059-016-0881-8"
            target="_blank">Conesa et al. 2016</a>).

            For RNA-seq QC analysis, QualiMap can be used to assess alignments produced by
            the first of these approaches. For input, it requires a GTF annotation file
            along with a reference genome, which can be used to reconstruct the exon
            structure of known transcripts. QualiMap uses this information to calculate the
            depth of coverage along the length of each annotated transcript. For a set of
            reads mapped to a transcript, the depth of coverage at a given base position is
            the number of high-quality reads that map to the transcript at that position
            (<a href="https://doi.org/10.1038/nrg3642" target="_blank">Sims et al. 2014</a>).

            QualiMap calculates coverage depth at every base position of each annotated
            transcript. To enable meaningful comparison between transcripts, base positions
            are rescaled to relative positions expressed as percentage distance along each
            transcript (*0%, 1%, &#8230;, 99%*). For the set of transcripts with at least
            one mapped read, QualiMap plots the cumulative mapped-read depth (y-axis) at
            each relative transcript position (x-axis). This plot shows the gene coverage
            profile across all mapped transcripts for each read dataset. It provides a
            visual way to assess positional biases, such as an accumulation of mapped reads
            at the 3&#8242; end of transcripts, which may indicate poor RNA quality in the
            original sample (<a href="https://doi.org/10.1186/s13059-016-0881-8"
            target="_blank">Conesa et al. 2016</a>).'''
            self.add_section(
                name='Gene Coverage Profile',
                anchor='qualimap-genome-fraction-coverage',
                description=
                'Mean distribution of coverage depth across the length of all mapped transcripts.',
                helptext=coverage_profile_helptext,
                plot=linegraph.plot(
                    self.qualimap_rnaseq_cov_hist, {
                        'id': 'qualimap_gene_coverage_profile',
                        'title':
                        'Qualimap RNAseq: Coverage Profile Along Genes (total)',
                        'ylab': 'Coverage',
                        'xlab': 'Transcript Position (%)',
                        'ymin': 0,
                        'xmin': 0,
                        'xmax': 100,
                        'tt_label': '<b>{point.x} bp</b>: {point.y:.0f}%',
                    }))

        # Return the number of reports we found
        return len(self.qualimap_rnaseq_genome_results.keys())
