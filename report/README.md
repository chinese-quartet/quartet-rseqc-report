---

![MultiReport for RNA-Seq](multireport-logo.png)

---

# MultiReport for RNA-Seq Pipeline

### Usage

To use this code, you need to install MultiQC and then your code. For example:

```bash
# Environment Configuration
conda create -n quartet-rseqc-report
conda activate quartet-rseqc-report
git clone https://github.com/chinese-quartet/quartet-rnaseq-report.git
cd ./quartet-rseqc-report/report
pip install -e .

# Usage
multiqc ./ (results direction)

# More Usage
multiqc -h
```

Use `python setup.py develop` if you're actively working on the code - then you don't need to rerun the installation every time you make an edit _(though you still do if you change anything in `setup.py`)_.

### Cautions

- post_alignment_qc/rnaseq_qc/ and post_alignment_qc/bam_qc/. All file names in both directories need to be in the form of samples_id.percent (i.e., the original form of qualimap output, without suffixes such as bam and rnase), and the IDs in both directories for the same sample should be identical
