"""Gene to OMIM Features V0.2

=== TCAG Summer Project 2018 ===
Yoonsik Park
BCB330Y1
University of Toronto

=== Module Description ===
Given a certain human entrez ID, this module will return OMIM features
This module also supports lookups with multiple human entrez IDs
"""

from typing import Dict, Tuple, List, Optional
import csv
import linecache

OMIM_FILENAME = "OMIMdiseasePhenotypeOnly_gene_info_GP20180213.tsv"
DELIMITER = "\t"
#Column numbers starting from zero
ENTREZ_COLUMN = 0
OMIM_DESCRIP_COLUMN = 4

class Gene_to_OMIM_Features:
    """This class loads the OMIM tsv file efficiently and allows for
    quick lookups of OMIM features

    === Attributes ===
    filename:
        The filename containing the MPO data. Most likely: 
        "OMIMdiseasePhenotypeOnly_gene_info_GP20180213.tsv"
    delimiter:
        The delimiter to be used. (i.e. comma or tab)

    === Private attributes ===
    _dict_genes_to_line:
        A dictionary to map entrez gene IDs to the corresponding
        line numbers of the file.
    """
    filename: str
    delimiter: str
    _dict_genes_to_line: Dict[int, int]

    def __init__(
            self,
            filename: str = OMIM_FILENAME,
            delimiter: str = DELIMITER
            ) -> None:
        """Opens the OMIM file, and creates data structures based on it."""
        self.filename = filename
        self.delimiter = delimiter
        self._dict_genes_to_line = {}
        file = open(filename, "r")
        file_reader = csv.reader(file, delimiter = delimiter)
        firstline = True
        line_number = 0
        for row in file_reader:
            line_number += 1
            if firstline:
                firstline = False
                continue
            entrez_id = int(row[ENTREZ_COLUMN])
            if entrez_id in self._dict_genes_to_line:
                print(__name__ + ": Error: multiple lines for the same entrez id!")
                print(__name__ + ": Row " + str(line_number) + ": " + str(row))
            else:
                self._dict_genes_to_line[entrez_id] = line_number
        file.close()

    def entrez_in_omim(self, entrez_id: int) -> bool:
        """Returns True if the entrez id is listed in the OMIM file
        Returns False otherwise
        """
        return entrez_id in self._dict_genes_to_line

    def entrez_to_omim_description(self, entrez_id: int) -> str:
        """Returns the description for the OMIM disease associated with the 
        entrez id. Raises a KeyError if the entrez id does not exist. You must 
        check if the entrez id is in OMIM before calling.
        """
        split_line = linecache.getline(self.filename,
                    self._dict_genes_to_line[entrez_id]).split(self.delimiter)
        return split_line[OMIM_DESCRIP_COLUMN] 

    def multi_entrez_to_num_diseases(self, entrez_ids: List[int]) -> int:
        """Returns the number of OMIM diseases that are associated with the 
        entrez ids. Returns 0 if no diseases are found.
        """
        omims_affected = 0
        for entrez_id in entrez_ids:
            if self.entrez_in_omim(entrez_id):
                omims_affected += 1
        return omims_affected

if __name__ == "__main__":
    omim = Gene_to_OMIM_Features()
    while True:
        try:
            input_string = input("Enter entrez id: ")
        except (EOFError, KeyboardInterrupt) as e:
            print("Exiting...")
            exit()
        try:
            entrez_id = int(input_string)
            print("Entrez ID in OMIM? "
                + str(omim.entrez_in_omim(entrez_id)))
        except ValueError:
            print("Invalid input")
        









