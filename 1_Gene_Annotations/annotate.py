from intervaltree import Interval, IntervalTree
import sys
import csv
import pickle
arguments = sys.argv[1:]
if len(arguments) != 2:
    print("annotate.py: requires")
    print("usage: python3 annotate.py [csvfile] [hg19_ncbi_refseq.txt.bin]")
    exit()

pickled_file = open(arguments[1], mode = 'rb')
list_chroms = pickle.load(pickled_file)

csvfile = open(arguments[0], newline = '')
outfile  = open(arguments[0].strip('.csv') + '_annotated.csv', "w+", newline = '')
readCSV = csv.DictReader(csvfile, delimiter = ',')
fieldnames = readCSV.fieldnames + ['genes_in_proximity', 'number_of_genes_in_proximity']
writer = csv.DictWriter(outfile, fieldnames)
writer.writeheader()

# chrX -> 22, chrY -> 23
chrom_dict = {'X': 22, 'Y': 23}
# chr1 -> 0, chr2 -> 1
for chrom_num in range(22):
    chrom_dict[str(chrom_num + 1)] = chrom_num

for row in readCSV:
    if 'chr' in row['chr']:
        chrom_num = row['chr'][3:]
    else:
        chrom_num = row['chr']
    if chrom_num in chrom_dict.keys():
        string_of_genes_affected = ""
        num_genes = 0
        for interval_ret in list_chroms[chrom_dict[chrom_num]].search(int(row['start']), int(row['end'])):
            string_of_genes_affected += interval_ret.data[0] + ";"
            num_genes += 1
        row['genes_in_proximity'] = string_of_genes_affected.strip(";")
        row['number_of_genes_in_proximity'] = str(num_genes)
        writer.writerow(row)
