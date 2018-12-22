"""Entrez to Gene Symbol Lookup V0.2

=== TCAG Summer Project 2018 ===
Yoonsik Park
BCB330Y1
University of Toronto

=== Module Description ===
Given either a certain human entrez ID or Gene Symbol, this module
will translate to the other gene identifier
"""

from typing import Dict, Tuple, List, Optional
import csv
import linecache

LOOKUP_FILENAME = "map.gsymbol.to.enzid.tsv"
DELIMITER = "\t"

class Entrez_Symbol_Lookup:
    """This class loads the entrez ID and gene symbol tsv file efficiently
    and allows for quick translation between either

    === Attributes ===
    filename:
        The filename containing the MPO data. Most likely: 
        "map.gsymbol.to.enzid.tsv"
    delimiter:
        The delimiter to be used. (i.e. comma or tab)
    dict_entrez_to_symbol:
        A dictionary to map entrez gene IDs to the corresponding
        symbol.
    dict_symbol_to_entrez:
        A dictionary to map gene symbols to the corresponding
        entrez gene ID.
    """
    filename: str
    delimiter: str
    dict_entrez_to_symbol: Dict[int, List[str]]
    dict_symbol_to_entrez: Dict[str, List[int]]

    def __init__(
            self,
            filename: str = LOOKUP_FILENAME,
            delimiter: str = DELIMITER
            ) -> None:
        """Opens the gene symbol/entrez file, and creates data structures based on it."""
        self.filename = filename
        self.delimiter = delimiter
        self.dict_entrez_to_symbol = {}
        self.dict_symbol_to_entrez = {}
        file = open(filename, "r")
        file_reader = csv.reader(file, delimiter = delimiter)
        firstline = True
        line_number = 0
        for row in file_reader:
            line_number += 1
            if firstline:
                firstline = False
                continue

            gene_symbol = row[0].strip().upper()
            entrez_id = int(row[1].strip())
            # ignore bad gene symbols
            if gene_symbol == '-':
                continue

            if entrez_id in self.dict_entrez_to_symbol:
                self.dict_entrez_to_symbol[entrez_id].append(gene_symbol)
            else:
                self.dict_entrez_to_symbol[entrez_id] = [gene_symbol]

            if gene_symbol in self.dict_symbol_to_entrez:
                self.dict_symbol_to_entrez[gene_symbol].append(entrez_id)
            else:
                self.dict_symbol_to_entrez[gene_symbol] = [entrez_id]
        counter = 0
        for value in self.dict_symbol_to_entrez.values():
            if len(value) != 1:
                counter += 1
        print(__name__ + ": Warning " + str(counter) + " symbol(s) have multiple entrez ids")
        counter =0
        for value in self.dict_entrez_to_symbol.values():
            if len(value) != 1:
                counter += 1
        print(__name__ + ": Warning " + str(counter)+ " entrez id(s) have multiple symbols")
        file.close()

    def entrez_id_to_symbols(self, entrez_id: int) -> List[str]:
        """Returns a list of symbols that are translations of the entrez id.
        Returns an empty list if the entrez id is not found.
        """
        if entrez_id not in self.dict_entrez_to_symbol:
            return []
        return(self.dict_entrez_to_symbol[entrez_id])

    def symbol_to_entrez_ids(self, symbol: str) -> List[int]:
        """Returns a list of entrez ids that are translations of the gene symbol
        Returns an empty list if the gene symbol is not found.
        """
        if symbol not in self.dict_symbol_to_entrez:
            return []
        return(self.dict_symbol_to_entrez[symbol])


if __name__ == "__main__":
    entrez_symbol = Entrez_Symbol_Lookup()
    while True:
        try:
            selection = input("Choose 'entrez' or 'symbol': ")
        except (EOFError, KeyboardInterrupt) as e:
            print("Exiting...")
            exit()
        try:
            my_id = input("Enter id: ")
        except (EOFError, KeyboardInterrupt) as e:
            print("Exiting...")
            exit()

        if selection[0].upper() == 'E':
            try:
                entrez_id = int(my_id)
                print("Symbol(s) found: "
                    + str(entrez_symbol.entrez_id_to_symbols(entrez_id)))
            except ValueError:
                print("Invalid input")
        else:
            try:
                symbol = my_id
                type(my_id)
                print("Entrez id(s) found: "
                    + str(entrez_symbol.symbol_to_entrez_ids(symbol)))
            except ValueError:
                print("Invalid input")











