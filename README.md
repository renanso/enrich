# enRich
Pipeline to design and test probes for target sequencing.

The first part of the pipeline is a script that extracts the flanking sequences of target SNPs and produces specific fragments. The second part performs sliding window to produce 120 bp probes. Finally BLAST is performed on the candidate probes to identify the potential binding sites and off targets. The third part takes the BLAST results and the sequences themselves to perform filtering for uniqueness, GC, Tm and thermodynamic calculations such as haipin and homodimer to select probes for multiplexing.

## Installation

```
install.packages("devtools")
devtools::install_github("renanso/enrich")
```
Key packages:
primer3
https://github.com/jensenlab/primer3

rBLAST
https://github.com/mhahsler/rBLAST

Biostrings
https://www.bioconductor.org/packages/release/bioc/html/Biostrings.html

dplyr
https://dplyr.tidyverse.org/

stringr
https://stringr.tidyverse.org/index.html

ggplot2
https://ggplot2.tidyverse.org/

## BLAST

For the package rBLAST to work, the BLAST software needs to be installed https://www.ncbi.nlm.nih.gov/books/NBK569861/

Add the BLAST executable to the PATH:
```
PATH=$PATH:~/path_to_folder/blast/ncbi-blast-2.XX.X+/bin
```
