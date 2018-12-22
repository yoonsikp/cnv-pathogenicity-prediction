import csv
import sys
arguments = sys.argv[1:]
if len(arguments) != 1:
    exit()
ofile  = open(arguments[0].split('.')[0] + '_fixed.csv', "w+", newline='')
writer = csv.writer(ofile)
file = open(arguments[0], 'r')
csvfile = csv.reader(file, delimiter=',')
firstline = True
for row in csvfile:
    if firstline:
        firstline = False
        writer.writerow(row)
        continue
    row[2] = 'chr' + row[2]
    writer.writerow(row)