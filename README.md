---

![MultiReport for RNA-Seq](multireport-logo.png)

---

# MultiReport for RNA-Seq Pipeline

### Usage

To use this code, you need to install MultiQC and then your code. For example:

```bash
# 环境配置
conda create -n multiqc multiqc
conda activate multiqc
git clone https://github.com/clinico-omics/chinese-quartet-rnaseq-report.git
cd rnaseq-report
pip install -r requirements.txt
pip install -e .

# 使用Multiqc命令
multiqc -h
```

Use `python setup.py develop` if you're actively working on the code - then you don't need to rerun the installation every time you make an edit _(though you still do if you change anything in `setup.py`)_.

### 注意事项

- post_alignment_qc/rnaseq_qc/和 post_alignment_qc/bam_qc/，这两个目录下所有文件名需要 samples_id.percent 形式，并且同一个样本两个目录下的 ID 要完全一致
