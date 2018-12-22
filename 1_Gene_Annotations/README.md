This folder contains the code used to convert the (chromosome, start, end) values to gene annotations.

### List of Files
* `hg19_ncbi_refseq.txt`

  Gene annotation database downloaded from NCBI RefSeq website. Using the hg19 genome assembly.

* `create_annotation_tree_entrez.py`

  This program assembles an IntervalTree for each chromosome with each gene's start and end location being an Interval. Each Interval contains the Entrez ID corresponding to the gene interval. This program will output a pickled (binary) file consisting of every IntervalTree. The resulting binary file is required for running the following annotation program.

* `hg19_ncbi_refseq.txt.bin`

  The pickled binary file containing all the gene annotation in Entrez ID format. This is explained above.

* `annotate.py`

  This program will perform the actual annotation on the CNV. The program requires a "chr", "start", and "end" column for the csv. This is the usage of the program:

   ```
   python3 annotate.py [input.csv] [hg19_ncbi_refseq.txt.bin]
   ```

  The program will then save a file called `[input]_annotated.csv` in the working directory
