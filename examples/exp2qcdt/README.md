> NOTE:
>
> 1. GENE_ID is the gene ID column name in the expression table, it must be there. And you need to make sure that the gene ID is unique and use the ENSEMBL gene ID which don't have version number. e.g. ENSG00000000003, not ENSG00000000003.14
>
> 2. The `fpkm.csv` and `count.csv` are the expression table, and the `metadata.csv` is the phenotype table.
>
> 3. The following columns are required in the phenotype table: `library`, `group`, `sample`, and the library ID must be identical to the expression table. The sample columns must be one of the following: `D5`, `D6`, `F7`, `M8`. The group column must be one of the following: `D5_1`, `D5_2`, `D5_3`, `D6_1`, `D6_2`, `D6_3`, `F7_1`, `F7_2`, `F7_3`, `M8_1`, `M8_2`, `M8_3`.