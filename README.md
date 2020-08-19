---

![MultiReport for RNA-Seq](multireport-logo.png)

---

# MultiReport for RNA-Seq Pipeline

### Usage

To use this code, you need to install MultiQC and then your code. For example:

```bash
# 环境配置
conda create -n multiqc multiqc
git clone https://github.com/clinico-omics/chinese-quartet-rnaseq-report.git
cd rnaseq-report
pip install -r requirements.txt
pip install -e .

# 使用Multiqc命令
multiqc -h
```

Use `python setup.py develop` if you're actively working on the code - then you don't need to rerun the installation every time you make an edit _(though you still do if you change anything in `setup.py`)_.
# quartet_rnaseq_qc_report
# chinese-quartet-rnaseq-report
