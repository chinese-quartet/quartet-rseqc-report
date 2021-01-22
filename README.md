# EXP2QCDT

Convert expression table to qc data table for quartet project.

## Installation

```R
devtools::install_github("clinico-omics/exp2qcdt")
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

## reference data

1.refqc_202011_forplot.rds
记录了 21 批次参考数据，6 种组合(D6/D5, F7/D5, M8/D5, F7/D6, M8/D6, M8/F7) relative correlation 值，以及 6 种组合参考数据集 logfc 之间的相关性，总共就是 21X6 行。第 4 列 corr_ref (Reference datasets based relative correlation) 就是参考数据集 logfc 之间的相关性，第 6 列 corr_FC (Relative correlation) 就是

> Reference datasets based relative correlation 解释
> 比如说：D6/D5，所有基因有一个 logfc，然后 D6/D5 与 21 个 batch D6/D5 logfc 做相关性，结果就是新的数据集 D6/D5 logfc 与参考数据之间的 logfc 相关性
> relative correlation 的解释
> 比如说：D6/D5，1. 三个 replicate 之间就有（D6_1/D5_1, D6_2/D5_1, D6_3/D5_1, D6_1/D5_2, D6_2/D5_2, D6_3/D5_2, D6_1/D5_3, D6_2/D5_3, D6_3/D5_3）9 种相对表达值；2. 9 种相对表达两两 correlation，就有了 36 个 correlation，然后取均值，就是每一个 batch，D6/D5 的 relative correlation

2. TableS2_ReferenceDatasets.csv
   > 6 种组合(D6/D5, F7/D5, M8/D5, F7/D6, M8/D6, M8/F7)，11797 个基因 logfc 的均值。

## Contributors

- [Jun Shang](https://github.com/stead99)
- [Jingcheng Yang](https://github.com/yjcyxky)
