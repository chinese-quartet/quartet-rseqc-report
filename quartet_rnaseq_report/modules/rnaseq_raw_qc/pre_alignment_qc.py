#!/usr/bin/env python
""" RnaSeqReport plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd
import numpy as np
import os
import zipfile
import re
import json
from collections import defaultdict
import plotly.express as px
import plotly.figure_factory as ff

from multiqc import config
from multiqc.plots import scatter, linegraph, table
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
            name='Pre Alignment QC',
            target="pre_alignment_qc",
            anchor='pre_alignment_qc',
            href='https://github.com/clinico-omics/quartet-rnaseq-report',
            info=" is an report module to show the quality of base.")

        self.pre_alignment_data_stats = list()
        self.pre_alignment_headers_stats = list()
        self.fastqc_data = dict()

        # Find and parse unzipped FastQC reports
        for f in self.find_log_files('rnaseq_raw_qc/data'):
            s_name = self.clean_s_name(os.path.basename(f['root']),
                                       os.path.dirname(f['root']))
            self.parse_fastqc_report(f['f'], s_name, f)

        # Find and parse zipped FastQC reports
        for f in self.find_log_files('rnaseq_raw_qc/zip', filecontents=False):
            s_name = f['fn']
            if s_name.endswith('_fastqc.zip'):
                s_name = s_name[:-11]
            # Skip if we already have this report - parsing zip files is slow..
            if s_name in self.fastqc_data.keys():
                log.debug("Skipping '{}' as already parsed '{}'".format(
                    f['fn'], s_name))
                continue
            try:
                fqc_zip = zipfile.ZipFile(os.path.join(f['root'], f['fn']))
            except Exception as e:
                log.warning("Couldn't read '{}' - Bad zip file".format(
                    f['fn']))
                log.debug("Bad zip file error:\n{}".format(e))
                continue
            # FastQC zip files should have just one directory inside, containing report
            d_name = fqc_zip.namelist()[0]
            try:
                with fqc_zip.open(os.path.join(d_name,
                                               'fastqc_data.txt')) as fh:
                    r_data = fh.read().decode('utf8')
                    self.parse_fastqc_report(r_data, s_name, f)
            except KeyError:
                log.warning(
                    "Error - can't find fastqc_raw_data.txt in {}".format(f))

        # Filter to strip out ignored sample names
        self.fastqc_data = self.ignore_samples(self.fastqc_data)

        if len(self.fastqc_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.fastqc_data)))
        # Write the summary stats to a file
        data = dict()
        for s_name in self.fastqc_data:
            data[s_name] = self.fastqc_data[s_name]['basic_statistics']
            data[s_name].update(self.fastqc_data[s_name]['statuses'])
        self.write_data_file(data, 'pre_alignment_fastqc')

        # Find and load any FastQ Screen reports
        self.fq_screen_data = dict()
        self.num_orgs = 0
        for f in self.find_log_files('rnaseq_raw_qc/fastq_screen',
                                     filehandles=True):
            parsed_data = self.parse_fqscreen(f)
            f['s_name'] = f['s_name'][:-7]
            if parsed_data is not None:
                if f['s_name'] in self.fq_screen_data:
                    log.debug(
                        "Duplicate sample name found! Overwriting: {}".format(
                            f['s_name']))
                self.add_data_source(f)
                self.fq_screen_data[f['s_name']] = parsed_data

        # Filter to strip out ignored sample names
        self.fq_screen_data = self.ignore_samples(self.fq_screen_data)

        if len(self.fq_screen_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.fq_screen_data)))

        # Check whether we have a consistent number of organisms across all samples
        num_orgs = set([len(orgs) for orgs in self.fq_screen_data.values()])
        # Write the total counts and percentages to files
        self.write_data_file(self.parse_csv(), 'pre_alignment_fastq_screen')

        # Add to self.css and self.js to be included in template
        self.css = {
            'assets/css/multiqc_fastqc.css':
            os.path.join(os.path.dirname(__file__), 'assets', 'css',
                         'multiqc_fastqc.css')
        }
        self.js = {
            'assets/js/multiqc_fastqc.js':
            os.path.join(os.path.dirname(__file__), 'assets', 'js',
                         'multiqc_fastqc.js')
        }

        # Colours to be used for plotting lines
        self.status_colours = {
            'pass': '#5cb85c',
            'warn': '#f0ad4e',
            'fail': '#d9534f',
            'default': '#999'
        }

        # Add the statuses to the intro for multiqc_fastqc.js JavaScript to pick up
        statuses = dict()
        for s_name in self.fastqc_data:
            for section, status in self.fastqc_data[s_name]['statuses'].items(
            ):
                try:
                    statuses[section][s_name] = status
                except KeyError:
                    statuses[section] = {s_name: status}
        self.intro += '<script type="application/json" class="fastqc_passfails">{}</script>'.format(
            json.dumps([self.anchor.replace('-', '_'), statuses]))

        self.intro += '<script type="text/javascript">load_fastqc_passfails();</script>'

        # Now add each section in order
        self.pre_aligment_stats()
        self.sequence_quality_plot()

    def parse_fastqc_report(self, file_contents, s_name=None, f=None):
        """ Takes contents from a fastq_data.txt file and parses out required
        statistics and data. Returns a dict with keys 'stats' and 'data'.
        Data is for plotting graphs, stats are for top table. """

        # Make the sample name from the input filename if we find it
        fn_search = re.search(r"Filename\s+(.+)", file_contents)
        if fn_search:
            s_name = self.clean_s_name(fn_search.group(1), f['root'])

        if s_name in self.fastqc_data.keys():
            log.debug(
                "Duplicate sample name found! Overwriting: {}".format(s_name))
        self.add_data_source(f, s_name)
        self.fastqc_data[s_name] = {'statuses': dict()}

        # Parse the report
        section = None
        s_headers = None
        self.dup_keys = []
        for l in file_contents.splitlines():
            if l == '>>END_MODULE':
                section = None
                s_headers = None
            elif l.startswith('>>'):
                (section, status) = l[2:].split("\t", 1)
                section = section.lower().replace(' ', '_')
                self.fastqc_data[s_name]['statuses'][section] = status
            elif section is not None:
                if l.startswith('#'):
                    s_headers = l[1:].split("\t")
                    # Special case: Total Deduplicated Percentage header line
                    if s_headers[0] == 'Total Deduplicated Percentage':
                        self.fastqc_data[s_name]['basic_statistics'].append({
                            'measure':
                            'total_deduplicated_percentage',
                            'value':
                            float(s_headers[1])
                        })
                    else:
                        # Special case: Rename dedup header in old versions of FastQC (v10)
                        if s_headers[1] == 'Relative count':
                            s_headers[1] = 'Percentage of total'
                        s_headers = [
                            s.lower().replace(' ', '_') for s in s_headers
                        ]
                        self.fastqc_data[s_name][section] = list()

                elif s_headers is not None:
                    s = l.split("\t")
                    row = dict()
                    for (i, v) in enumerate(s):
                        v.replace('NaN', '0')
                        try:
                            v = float(v)
                        except ValueError:
                            pass
                        row[s_headers[i]] = v
                    self.fastqc_data[s_name][section].append(row)
                    # Special case - need to remember order of duplication keys
                    if section == 'sequence_duplication_levels':
                        try:
                            self.dup_keys.append(float(s[0]))
                        except ValueError:
                            self.dup_keys.append(s[0])

        # Tidy up the Basic Stats
        self.fastqc_data[s_name]['basic_statistics'] = {
            d['measure']: d['value']
            for d in self.fastqc_data[s_name]['basic_statistics']
        }

        # Calculate the average sequence length (Basic Statistics gives a range)
        length_bp = 0
        total_count = 0
        for d in self.fastqc_data[s_name].get('sequence_length_distribution',
                                              {}):
            length_bp += d['count'] * self.avg_bp_from_range(d['length'])
            total_count += d['count']
        if total_count > 0:
            self.fastqc_data[s_name]['basic_statistics'][
                'avg_sequence_length'] = length_bp / total_count

    def pre_aligment_stats(self):
        """ Add some single-number stats to the basic statistics
        table at the top of the report """

        # Prep the data
        data = dict()
        for s_name in self.fastqc_data:
            bs = self.fastqc_data[s_name]['basic_statistics']
            try:
                # FastQC reports with 0 reads will trigger a KeyError here
                data[s_name] = {
                    'total_sequences': bs['Total Sequences'],
                    'percent_gc': bs['%GC']
                }
            except KeyError:
                log.warning("Sample had zero reads: '{}'".format(s_name))
                data[s_name] = {'total_sequences': 0, 'percent_gc': 0}
        ts_mean = []
        gc_mean = []
        for k in data.keys():
            ts_mean.append(data[k]['total_sequences'])
            gc_mean.append(data[k]['percent_gc'])
        data['Batch average value'] = {
            'total_sequences': sum(ts_mean) / len(ts_mean),
            'percent_gc': sum(gc_mean) / len(gc_mean)
        }
        data['Historical value'] = {
            'total_sequences': 59530000,
            'percent_gc': 48.82
        }

        headers = OrderedDict()
        headers['total_sequences'] = {
            'title': '{} Seqs'.format(config.read_count_prefix),
            'description':
            'Total Sequences ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        headers['percent_gc'] = {
            'title': '% GC',
            'description': 'Average % GC Content',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Set1',
            'format': '{:,.1f}'
        }
        headers['Human percentage'] = {
            'title': '% Human',
            'description': 'Average % Human',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Set1',
            'format': '{:,.0f}'
        }
        self.pre_alignment_data_stats.append(data)
        self.pre_alignment_data_stats.append(self.parse_csv())
        self.pre_alignment_headers_stats.append(headers)
        self.add_section(name='Pre alignment stats',
                         anchor='pre_alignment_stats',
                         description='The summary of fastqc.',
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
                         plot=table.plot(self.pre_alignment_data_stats,
                                         self.pre_alignment_headers_stats))

    def sequence_quality_plot(self):
        """ Create the HTML for the phred quality score plot """

        data = dict()
        for s_name in self.fastqc_data:
            try:
                data[s_name] = {
                    self.avg_bp_from_range(d['base']): d['mean']
                    for d in self.fastqc_data[s_name]
                    ['per_base_sequence_quality']
                }
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('sequence_quality not found in FastQC reports')
            return None

        pconfig = {
            'id':
            'fastqc_per_base_sequence_quality_plot',
            'title':
            'FastQC: Mean Quality Scores',
            'ylab':
            'Phred Score',
            'xlab':
            'Position (bp)',
            'ymin':
            0,
            'xDecimals':
            False,
            'tt_label':
            '<b>Base {point.x}</b>: {point.y:.2f}',
            'colors':
            self.get_status_cols('per_base_sequence_quality'),
            'yPlotBands': [
                {
                    'from': 28,
                    'to': 100,
                    'color': '#c3e6c3'
                },
                {
                    'from': 20,
                    'to': 28,
                    'color': '#e6dcc3'
                },
                {
                    'from': 0,
                    'to': 20,
                    'color': '#e6c3c3'
                },
            ]
        }
        self.add_section(
            name='Sequence Quality Histograms',
            anchor='fastqc_per_base_sequence_quality',
            description=
            'The mean quality value across each base position in the read.',
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
            plot=linegraph.plot(data, pconfig))

    def avg_bp_from_range(self, bp):
        """" Helper function - FastQC often gives base pair ranges (eg. 10-15)
        which are not helpful when plotting. This returns the average from such
        ranges as an int, which is helpful. If not a range, just returns the int. """

        try:
            if '-' in bp:
                maxlen = float(bp.split("-", 1)[1])
                minlen = float(bp.split("-", 1)[0])
                bp = ((maxlen - minlen) / 2) + minlen
        except TypeError:
            pass
        return (int(bp))

    def get_status_cols(self, section):
        """ Helper function - returns a list of colours according to the FastQC
        status of this module for each sample. """

        colours = dict()
        for s_name in self.fastqc_data:
            status = self.fastqc_data[s_name]['statuses'].get(
                section, 'default')
            colours[s_name] = self.status_colours[status]
        return colours

    def parse_fqscreen(self, f):
        """ Parse the FastQ Screen output into a 3D dict """
        parsed_data = OrderedDict()
        nohits_pct = None
        headers = None
        bs_headers = None
        for l in f['f']:
            # Skip comment lines
            if l.startswith('#'):
                continue
            if l.startswith('%Hit_no_genomes:') or l.startswith(
                    '%Hit_no_libraries:'):
                nohits_pct = float(l.split(':', 1)[1])
                parsed_data['No hits'] = {
                    'percentages': {
                        'one_hit_one_library': nohits_pct
                    }
                }
            else:
                s = l.strip().split("\t")

                # Regular FastQ Screen table section
                if len(s) == 12:
                    if headers is None:
                        headers = s
                    else:
                        # Can't use #Reads in subset as varies. #Reads_processed should be same for all orgs in a sample
                        parsed_data['total_reads'] = int(s[1])
                        # Loop through all columns
                        parsed_data[s[0]] = {'percentages': {}, 'counts': {}}
                        for idx, h in enumerate(headers):
                            if idx == 0:
                                continue
                            dtype = 'percentages' if h.startswith(
                                '%') else 'counts'
                            field = h.replace('%', '').replace(
                                '#',
                                '').replace('genomes', 'libraries').replace(
                                    'genome', 'library').lower()
                            parsed_data[s[0]][dtype][field] = float(s[idx])

                # Optional Bisulfite table section
                elif len(s) == 9:
                    if bs_headers is None:
                        bs_headers = s
                    else:
                        # Loop through all columns
                        parsed_data[s[0]]['bisulfite_percentages'] = {}
                        parsed_data[s[0]]['bisulfite_counts'] = {}
                        for idx, h in enumerate(bs_headers):
                            if idx == 0:
                                continue
                            dtype = 'bisulfite_percentages' if h.startswith(
                                '%') else 'bisulfite_counts'
                            field = h.replace('%', '').replace('#', '').lower()
                            parsed_data[s[0]][dtype][field] = float(s[idx])

        if len(parsed_data) == 0:
            return None

        # Calculate no hits counts
        if parsed_data.get(
                'total_reads') is not None and nohits_pct is not None:
            parsed_data['No hits']['counts'] = {
                'one_hit_one_library':
                int((nohits_pct / 100.0) * float(parsed_data['total_reads']))
            }
        else:
            log.warning(
                "Couldn't find number of reads with no hits for '{}'".format(
                    f['s_name']))

        self.num_orgs = max(len(parsed_data), self.num_orgs)
        return parsed_data

    def parse_csv(self):
        totals = dict()
        for s in sorted(self.fq_screen_data.keys()):
            totals[s] = dict()
            for org in self.fq_screen_data[s]:
                if org == 'total_reads':
                    continue
                try:
                    k = "{} percentage".format(org)
                    totals[s][k] = self.fq_screen_data[s][org]['percentages'][
                        'one_hit_one_library']
                    totals[s][k] += self.fq_screen_data[s][org][
                        'percentages'].get('multiple_hits_one_library', 0)
                    totals[s][k] += self.fq_screen_data[s][org][
                        'percentages'].get('one_hit_multiple_libraries', 0)
                    totals[s][k] += self.fq_screen_data[s][org][
                        'percentages'].get('multiple_hits_multiple_libraries',
                                           0)
                except KeyError:
                    pass
        human_mean = []
        ercc_mean = []
        ecoli_mean = []
        adapter_mean = []
        vector_mean = []
        rrna_mean = []
        virus_mean = []
        yeast_mean = []
        mitoch_mean = []
        no_hits_mean = []
        for k in totals.keys():
            human_mean.append(totals[k]['Human percentage'])
            ercc_mean.append(totals[k]['ERCC percentage'])
            ecoli_mean.append(totals[k]['EColi percentage'])
            adapter_mean.append(totals[k]['Adapter percentage'])
            vector_mean.append(totals[k]['Vector percentage'])
            rrna_mean.append(totals[k]['rRNA percentage'])
            virus_mean.append(totals[k]['Virus percentage'])
            yeast_mean.append(totals[k]['Yeast percentage'])
            mitoch_mean.append(totals[k]['Mitoch percentage'])
            no_hits_mean.append(totals[k]['No hits percentage'])
        totals['Batch average value'] = {
            'Human percentage': sum(human_mean) / len(human_mean),
            'ERCC percentage': sum(ercc_mean) / len(ercc_mean),
            'EColi percentage': sum(ecoli_mean) / len(ecoli_mean),
            'Adapter percentage': sum(adapter_mean) / len(adapter_mean),
            'Vector percentage': sum(vector_mean) / len(vector_mean),
            'rRNA percentage': sum(rrna_mean) / len(rrna_mean),
            'Virus percentage': sum(virus_mean) / len(virus_mean),
            'Yeast percentage': sum(yeast_mean) / len(yeast_mean),
            'Mitoch percentage': sum(mitoch_mean) / len(mitoch_mean),
            'No hits percentage': sum(no_hits_mean) / len(no_hits_mean),
        }
        totals['Historical value'] = {
            'Human percentage': '98.34±14.33',
            'ERCC percentage': '0.15±0.02',
            'EColi percentage': '0.003±0.004',
            'Adapter percentage': '0.003±0.001',
            'Vector percentage': '0.19±0.03',
            'rRNA percentage': '1.97±0.29',
            'Virus percentage': '0.64±0.09',
            'Yeast percentage': '0.34±0.05',
            'Mitoch percentage': '2.70±0.39',
            'No hits percentage': '0.99±0.14'
        }
        return totals
