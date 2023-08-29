# EXP2QCDT

Convert expression table to qc data table for quartet project.

## Installation

```R
devtools::install_github("chinese-quartet/quartet-rseqc-report", subdir = "exp2qcdt")
```

## Usage

```R
library(exp2qcdt)
exp2qcdt("~/Downloads/exp2qcdt/test/fpkm_table.txt", "~/Downloads/exp2qcdt/test/counts_table.txt", "~/Downloads/exp2qcdt/test/phenotype.txt", "~/Downloads/exp2qcdt/test/")
```

## phenotype table example

| library                                         | group | sample |
| :---------------------------------------------- | ----- | ------ |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_D5_1_20200618 | D5_1  | D5     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_D5_2_20200618 | D5_2  | D5     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_D5_3_20200618 | D5_3  | D5     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_D6_1_20200618 | D6_1  | D6     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_D6_2_20200618 | D6_2  | D6     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_D6_3_20200618 | D6_3  | D6     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_F7_1_20200618 | F7_1  | F7     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_F7_2_20200618 | F7_2  | F7     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_F7_3_20200618 | F7_3  | F7     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_M8_1_20200618 | M8_1  | M8     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_M8_2_20200618 | M8_2  | M8     |
| Quartet_RNA_BGI_BGI2000_PolyA_BGI_M8_3_20200618 | M8_3  | M8     |

The fpkm and counts library IDs need to be identical to phenotype

## reference data

> ref_data_fc_value.csv
> The log2fc of D5/D6, F7/D6, M8/D6 were calculated for the filtered genes based on the reference dataset.

> ref_data_qc_value.csv
> SNR, RC, Total Score value of 21 historical datasets.

## Contributors

- [Jun Shang](https://github.com/stead99)
- [Jingcheng Yang](https://github.com/yjcyxky)
