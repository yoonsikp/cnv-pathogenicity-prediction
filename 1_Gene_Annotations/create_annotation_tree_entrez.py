from intervaltree import Interval, IntervalTree
import csv
from Bio import Entrez
import pickle

Entrez.email = 'yoonsik.park@mail.utoronto.ca'
filename = 'hg19_ncbi_refseq.txt'
pickled_file = open(filename + ".bin", 'wb+')
csvfile = open(filename, newline = '')
readCSV = csv.reader(csvfile, delimiter = '\t')

# [{32334: (start, end), ... }, {22334: (start, end), ... }, ... every chromosome ... ]
all_genes = []

for _ in range(24):
    all_genes.append({})

# chrX -> chromosome 23, chrY -> chromosome 24
chrom_dict = {'X': 22, 'Y': 23}
# chr1 -> 0, chr2 -> 1
for chrom_num in range(22):
    chrom_dict[str(chrom_num + 1)] = chrom_num 

firstline = True
counter = 0
for row in readCSV:
    counter += 1
    print(counter)
    # skip first line
    if firstline:    
        firstline = False
        continue
    if counter == 15:
        break
    # check if chromosome is valid and if gene is valid
    chrom_num = row[2][3:]
    if chrom_num in chrom_dict.keys() and row[13].strip() != "none":
        # convert accession number to entrez ID
        accession_number = row[1].strip().split(".")[0]
        esearch_result = Entrez.esearch(db = "gene",
                                        term = (accession_number + "[ACCN]"),
                                        retmod = "xml")
        parsed_result = Entrez.read(esearch_result)
        if len(parsed_result['IdList']) != 1:
            print("entrez id couldn't be found." + str(parsed_result) + row[1])
            continue
        else:
            entrez_id = parsed_result['IdList'][0]

        # add gene interval to the right all_genes chromosome
        if entrez_id in all_genes[chrom_dict[chrom_num]].keys():
            gene_tuple = all_genes[chrom_dict[chrom_num]][entrez_id]
            # get widest interval
            all_genes[chrom_dict[chrom_num]][entrez_id] = min(gene_tuple[0], int(row[4])),\
                                                          max(gene_tuple[1], int(row[5])) 
        else:
            all_genes[chrom_dict[chrom_num]][entrez_id] = int(row[4]), int(row[5])

# convert to efficient tree structure
list_chroms = []
for i in range(24):
    list_chroms.append(IntervalTree())
    for key,value in all_genes[i].items():
        list_chroms[i][value[0]:value[1]] = (key, None)

print(list_chroms)
pickle.dump(list_chroms, pickled_file)
pickled_file.close()

