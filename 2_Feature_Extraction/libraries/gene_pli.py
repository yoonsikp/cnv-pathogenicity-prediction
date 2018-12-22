"""Gene to pLI Features V0.2

=== TCAG Summer Project 2018 ===
Yoonsik Park
BCB330Y1
University of Toronto

=== Module Description ===
Given a gene symbol, this module will return the pLI
This module is extendable to return other columns from the pLI file
"""

from typing import Dict, Tuple, List, Optional, Union
import csv
import linecache

PLI_FILENAME = "pLI_EXac_broad_institure_2016_03.txt"
DELIMITER = "\t"
#Column numbers starting from zero
GENE_SYMBOL_COLUMN = 1
PLI_COLUMN = 19

class Gene_to_pLI_Features:
    """This class loads the pLI file efficiently and allows for
    quick lookups of pLI features based on gene symbol.

    === Attributes ===
    filename:
        The filename containing the pLI data. Most likely: 
        "pLI_EXac_broad_institure_2016_03.txt"
    delimiter:
        The delimiter to be used. (i.e. comma or tab)

    === Private attributes ===
    _dict_gene_symbol_to_line:
        A dictionary to map gene symbols to the corresponding
        line numbers of the file.
    """
    filename: str
    delimiter: str
    _dict_gene_symbol_to_line: Dict[str, int]

    def __init__(
            self,
            filename: str = PLI_FILENAME,
            delimiter: str = DELIMITER
            ) -> None:
        """Opens the pLI file, and creates data structures based on it."""
        self.filename = filename
        self.delimiter = delimiter
        self._dict_gene_symbol_to_line = {}
        file = open(filename, "r")
        file_reader = csv.reader(file, delimiter = delimiter)
        firstline = True
        line_number = 0
        for row in file_reader:
            line_number += 1
            if firstline:
                firstline = False
                continue
            symbol = row[GENE_SYMBOL_COLUMN].upper()
            if symbol in self._dict_gene_symbol_to_line:

                print(__name__ + ": Warning: multiple lines for the same gene symbol!")
                print(__name__ + ": Row " + str(line_number) + ": " + str(row))
            else:
                self._dict_gene_symbol_to_line[symbol] = line_number
        file.close()

    def gene_symbol_to_pLI(self, symbol: str) -> Union[float, int]:
        """Returns the floating point pLI of the gene symbol.
        Returns -1 if the gene symbol is not found
        """
        if symbol not in self._dict_gene_symbol_to_line:
            return -1
        return float(linecache.getline(self.filename, self._dict_gene_symbol_to_line[symbol]).split(self.delimiter)[PLI_COLUMN])

        # for any other column:
        # if symbol not in self._dict_gene_symbol_to_line:
        #     return -1
        # return float(linecache.getline(self.filename, self._dict_gene_symbol_to_line[symbol]).split(self.delimiter)[***REPLACE_HERE***])



if __name__ == "__main__":
    pli = Gene_to_pLI_Features()
    while True:
        try:
            symbol = input("Enter gene symbol: ").upper()
        except (EOFError, KeyboardInterrupt) as e:
            print("Exiting...")
            exit()
        pli_ret = pli.gene_symbol_to_pLI(symbol)
        if pli_ret != -1:
            print("pLI: " + str(pli_ret))
        else:
            print("Gene symbol not found!")








