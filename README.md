# ATG_CRISPResso2_To_RIMA2

## Background
Non-Homologous End Joining (NHEJ) and Micohomology Mediated End Joining (MMEJ) are two highly conserved DNA rapir pathways in eukaryotic cells. 

## Requirements

```
crispresso2=2.2.14
fastx_toolkit=0.0.14
flash=1.2.11
```

## Step-by-step guideline (Mac and Linux)

### Install Miniconda
Conda is an open-source package and environment management system that allows user to conviniently create new environments and install packages. Follow the official Conda installation guideline here: https://docs.anaconda.com/free/miniconda/miniconda-install/ 

### Download this repository
```
git clone https://github.com/ghahfarokhi/atg_crispresso2_to_rima2
cd atg_crispresso2_to_rima2
```

### Create atg_crispresso2_to_rima2 environment
```
conda env create -f environment.yml
conda activate atg_crispresso2_to_rima2
```

Test that the installation of CRISPResso2 and fastx_toolkit has been successful: 
* `CRISPResso -h` prints the help page for CRISPResso2.
* `fastx_barcode_splitter.pl -h` prints the help page for fastx_barcode_splitter.

For detailed description of CRISPResso functionalities and alternative ways of installation refer to its official github page: https://github.com/pinellolab/CRISPResso2

### Run CRISPResso
Run CRISPResso with a desired set of paramaters (`-a` and `-g` are required parameters for RIMA analysis). For this guideline, the following codes uses the test fastq data that are provided withinn this repository with the following amplicon sequence "CGAGTCTAGAGGGCCCGTTTAAACCCGCTGGGCCATGGGCTATGAGTACAGGTCATGTACGGCCTCATAGTGGTACAGTAGTGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTA" and guide sequence "CTATGAGTACAGGTCATGTA" (also provided in a fasta format in `test_data/01_fastq_examples/reference_and_guide_sequence.fasta`). 
```
CRISPResso -r1 01_fastq_examples/example_DMSO.fastq -a CGAGTCTAGAGGGCCCGTTTAAACCCGCTGGGCCATGGGCTATGAGTACAGGTCATGTACGGCCTCATAGTGGTACAGTAGTGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTA -g CTATGAGTACAGGTCATGTA -o 02_crispresso_examples/

CRISPResso -r1 01_fastq_examples/example_NHEJi.fastq -a CGAGTCTAGAGGGCCCGTTTAAACCCGCTGGGCCATGGGCTATGAGTACAGGTCATGTACGGCCTCATAGTGGTACAGTAGTGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTA -g CTATGAGTACAGGTCATGTA -o 02_crispresso_examples/
```

The test fastq files provided in `test_data/01_fastq_examples/` as example are simulated for two samples transfected with SpyCas9 and a sgRNA, and treated with either DMSO or an NHEJ-inhibitory compound (NHEJi). 

See CRISPResso output for DMSO and NHEJi samples:
* [CRISPResso_on_example_DMSO.html](https://ghahfarokhi.github.io/ATG_CRISPResso2_To_RIMA2/test_data/02_crispresso_examples/02_crispresso_examples/CRISPResso_on_example_DMSO.html)
* [CRISPResso_on_example_NHEJi.html](https://ghahfarokhi.github.io/ATG_CRISPResso2_To_RIMA2/test_data/02_crispresso_examples/02_crispresso_examples/CRISPResso_on_example_NHEJi.html)

### Run `ATG_CRISPResso2_to_RIMA2.py` script
``` 
python ATG_CRISPResso2_to_RIMA2.py --crispresso_folder 02_crispresso_examples/ --out 03_rima_examples/
```

Output:

```
example_DMSO
	Processing allele frequency table: 02_crispresso_examples/CRISPResso_on_example_DMSO/Alleles_frequency_table.zip
	Variant file for example_DMSO saved as 03_rima_examples/RIMA_example_DMSO-variants.tsv
example_NHEJi
	Processing allele frequency table: 02_crispresso_examples/CRISPResso_on_example_NHEJi/Alleles_frequency_table.zip
	Variant file for example_NHEJi saved as 03_rima_examples/RIMA_example_NHEJi-variants.tsv

Tab-delimited experiment_sheet file saved as 03_rima_examples/experiment_sheet.tsv
Done!
```

The following three files should have been generated in the output folder upon successful completion of the python script run:

`ls 03_rima_examples/`

Output:
``` 
RIMA_example_DMSO-variants.tsv
RIMA_example_NHEJi-variants.tsv
experiment_sheet.tsv
```

`column -t 03_rima_examples/RIMA_example_DMSO-variants.tsv`

Output:
```
VariantNo  Position  Type  Length     Ref  Alt           Count
0          2         57    Insertion  1    -             T      35
1          3         50    Deletion   11   AGGTCATGTAC   -      25
2          4         54    Deletion   5    CATGT         -      16
3          5         57    Deletion   12   GTACGGCCTCAT  -      10
4          6         53    Deletion   10   TCATGTACGG    -      9
5          7         55    Deletion   2    AT            -      8
6          8         57    Deletion   1    G             -      6
7          9         56    Deletion   1    T             -      3
8          10        56    Deletion   4    TGTA          -      2
```

`column -t 03_rima_examples/experiment_sheet.tsv`

Output:
```
file_address  file_name                                    amplicon_seq                     guide_seq                                                                                                                                                                                                   mapped_reads          wt_count  total_variants  cut_pos
0             C:\RIMA\Raw\RIMA_example_DMSO-variants.tsv   RIMA_example_DMSO-variants.tsv   CGAGTCTAGAGGGCCCGTTTAAACCCGCTGGGCCATGGGCTATGAGTACAGGTCATGTACGGCCTCATAGTGGTACAGTAGTGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTA  CTATGAGTACAGGTCATGTA  180       66              9        56
1             C:\RIMA\Raw\RIMA_example_NHEJi-variants.tsv  RIMA_example_NHEJi-variants.tsv  CGAGTCTAGAGGGCCCGTTTAAACCCGCTGGGCCATGGGCTATGAGTACAGGTCATGTACGGCCTCATAGTGGTACAGTAGTGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTA  CTATGAGTACAGGTCATGTA  180       77              6        56

```

### Transfer RIMA files to a Windows computer
* Copy and paste `experiment_sheet.tsv` table to `Experiment` worksheet in **RIMA2_v20240420.xlsm** file.
* Transfer RIMA variant tables to **"C:\RIMA\Raw\"** folder or adjust the file addresses accordign to your local disk address. *IMPORTANT NOTE*: file addresses should be alphanumerical without spaces.

### RIMA vs CRISPResso
![Snapshots of RIAM and CRISPResso outputs](./test_data/CRISPResso_vs_RIMA_output.png)