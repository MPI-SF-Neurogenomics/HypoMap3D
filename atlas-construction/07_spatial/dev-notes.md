# Spatial Development Notes

## Dataset Notes

### Yao | Zeng | MERFISH | ABC | SP-01

A few gene names were manually mapped, due to Gene Symbol changes, to match reference genes.

### Zhang | Zhuang | MERFISH | ABC | SP2-02

A few gene names were manually mapped, due to Gene Symbol changes, to match reference genes.

### Bauder | OpenST | MPI | SP2-03

Gene names were derived from Ensembl Version 107 - they were mapped to 110 to match single cell.

A large number of CCF coordinates fall out of bounds (OOB). Why?

### Langlieb | Macosko | Slide-seq | BrainCellData | SP2-04

Gene names were derived from Ensembl Version 103 - they were mapped to 110 to match single cell. 8 genes were lost in this process.

Several pucks had slides re-registered. For further details see 3d_src_Langlieb.ipynb.

#### Harmonization Notes

Several metadata were viable as the `Sample` column including `orig.ident`, `puck_num`, `Puck_Depth` and `PuckID`. `PuckID` was chosen as the best candidate as the most descriptive among them.

## Module Notes

### CCF Setup

SFO is defined as hypothalamic according to ABC taxonomy. This was manually adjusted to be excluded as non-HY. This comes with the significant advantage of reducing CCF Y maximum values.
