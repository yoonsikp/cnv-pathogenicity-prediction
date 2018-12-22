import csv
file = open("all_labels.csv", 'r')
wfile = open("all_labels.tsv", 'w')
csvreader = csv.reader(file)
csvwriter = csv.writer(wfile, delimiter='\t')
for row in csvreader:
    csvwriter.writerow(row)
