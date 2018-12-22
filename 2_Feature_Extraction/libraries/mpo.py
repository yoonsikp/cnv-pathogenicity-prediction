"""Gene to MPO Features V0.2

=== TCAG Summer Project 2018 ===
Yoonsik Park
BCB330Y1
University of Toronto

=== Module Description ===
Given a certain human entrez ID, this module will return MPO features
This module also supports lookups with multiple human entrez IDs
"""

from typing import Dict, Tuple, List, Optional
import csv
import linecache

MPO_FILENAME = "MPO_topPh_GP20180216.tsv"
DELIMITER = "\t"
#Column numbers starting from zero
ENTREZ_COLUMN = 0
PHENOTYPE_COLUMN = 4

class Gene_to_MPO_Features:
    """This class loads the MPO tsv file efficiently and allows for
    quick lookups of MPO features

    === Attributes ===
    filename:
        The filename containing the MPO data. Most likely: 
        "MPO_topPh_GP20180216.tsv"
    delimiter:
        The delimiter to be used. (i.e. comma or tab)

    === Private attributes ===
    _dict_genes_to_line:
        A dictionary to map entrez gene IDs to the corresponding
        line numbers of the file.
    _list_phenotypes:
        A sorted list of all phenotype strings
    """
    filename: str
    delimiter: str
    _dict_genes_to_line: Dict[int, List[int]]
    _list_phenotypes: List[str]

    def __init__(
            self,
            filename: str = MPO_FILENAME,
            delimiter: str = DELIMITER
            ) -> None:
        """Opens the mpo file, and creates data structures based on it."""
        self.filename = filename
        self.delimiter = delimiter
        self._dict_genes_to_line = {}
        self._list_phenotypes = []
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
                self._dict_genes_to_line[entrez_id].append(line_number)
            else:
                self._dict_genes_to_line[entrez_id] = [line_number]
            phenotype = row[PHENOTYPE_COLUMN]
            if phenotype not in self._list_phenotypes:
                self._list_phenotypes.append(phenotype)
        self._list_phenotypes.sort()
        file.close()

    def get_all_phenotypes(self) -> List[str]:
        """Returns the list of phenotypes found in the file."""
        return self._list_phenotypes

    def entrez_to_num_phenotypes(self, entrez_id: int) -> int:
        """Returns the number of phenotypes that are associated with the entrez id.
        This is equivalent to the number of lines in the file with the entrez id.
        Returns 0 if no phenotypes are found.
        """
        if entrez_id not in self._dict_genes_to_line:
            return 0
        return(len(self._dict_genes_to_line[entrez_id]))

    def entrez_to_num_uniq_phenotypes(self, entrez_id: int) -> int:
        """Returns the number of unique phenotypes that are associated with the entrez id.
        This function exists because there are repeated phenotypes for a single entrez id.
        The maximum returned value is len(self.get_all_phenotypes())
        Returns 0 if no phenotypes are found. 
        """
        if entrez_id not in self._dict_genes_to_line:
            return 0

        phenotypes_affected = 0
        for value in self.entrez_to_phenotypes(entrez_id).values():
            if value > 0:
                phenotypes_affected += 1
        return phenotypes_affected

    def entrez_to_phenotypes(self, entrez_id: int) -> Dict[str, int]:
        """Returns a dict of all phenotypes mapping to the number of times the
        entrez id was associated. Returns a dict with values of all zero
        if no phenotypes are found.
        """
        dict_phenotype_to_hits = {}

        for phenotype in self.get_all_phenotypes():
            dict_phenotype_to_hits[phenotype] = 0

        if entrez_id in self._dict_genes_to_line:
            for line_num in self._dict_genes_to_line[entrez_id]:
                phenotype = linecache.getline(
                    self.filename, line_num).split(self.delimiter)[PHENOTYPE_COLUMN]
                dict_phenotype_to_hits[phenotype] += 1

        return dict_phenotype_to_hits

    def thresh_entrez_to_phenotypes(self, entrez_id: int) -> Dict[str, int]:
        """Returns a dict of all phenotypes mapping to either a 1 or 0 depending if
        entrez id was associated or not respectively. Returns a dict with values of all zero
        if no phenotypes are found.
        """
        dict_phenotype_to_hits = self.entrez_to_phenotypes(entrez_id)

        for key, value in dict_phenotype_to_hits.items():
            if value > 1:
                dict_phenotype_to_hits[key] = 1

        return dict_phenotype_to_hits

    def multi_entrez_to_num_phenotypes(self, entrez_ids: List[int]) -> int:
        """Returns the number of phenotypes that are associated with the entrez
        ids. This is equivalent to the number of lines in the file that match any
        of the entrez ids.
        Returns 0 if no phenotypes are found.
        """
        phenotypes_affected = 0
        for entrez_id in entrez_ids:
            phenotypes_affected += self.entrez_to_num_phenotypes(entrez_id)
        return phenotypes_affected

    def multi_entrez_to_num_phenotypes_using_thresh(self, entrez_ids: List[int]) -> int:
        """Returns the number of phenotypes that are associated with the entrez
        ids. Uses thresh for each entrez id.
        Returns 0 if no phenotypes are found.
        """
        phenotypes_affected = 0
        for entrez_id in entrez_ids:
            phenotypes_affected += self.entrez_to_num_uniq_phenotypes(entrez_id)
        return phenotypes_affected

    def multi_entrez_to_num_uniq_phenotypes(self, entrez_ids: List[int]) -> int:
        """Returns the number of unique phenotypes that are associated with the entrez
        ids. This function exists because there are repeated phenotypes for a single entrez id.
        The maximum returned value is len(self.get_all_phenotypes()).
        Returns 0 if no phenotypes are found.
        """

        phenotypes_affected = 0
        for value in self.multi_entrez_to_phenotypes(entrez_ids).values():
            if value > 0:
                phenotypes_affected += 1
        return phenotypes_affected

    def multi_entrez_to_phenotypes(
            self,
            entrez_ids: List[int]
            ) -> Dict[str, int]:
        """Returns a dict of the phenotypes and number of hits that 
        are associated with the entrez ids.
        """
        dict_phenotype_to_hits = {}
        for phenotype in self.get_all_phenotypes():
            dict_phenotype_to_hits[phenotype] = 0
        for entrez_id in entrez_ids:
            for key, value in self.entrez_to_phenotypes(entrez_id).items():
                dict_phenotype_to_hits[key] += value
        return dict_phenotype_to_hits

    def thresh_multi_entrez_to_phenotypes(
            self,
            entrez_ids: List[int]
            ) -> Dict[str, int]:
        """Returns a dict of all phenotypes mapping to either a 1 or 0 depending if any of the
        entrez ids were associated or not, respectively. 
        """
        dict_phenotype_to_hits = {}
        for phenotype in self.get_all_phenotypes():
            dict_phenotype_to_hits[phenotype] = 0
        for key, value in self.multi_entrez_to_phenotypes(entrez_ids).items():
            if value > 0:
                dict_phenotype_to_hits[key] = 1
        return dict_phenotype_to_hits

    def multi_entrez_to_phenotypes_using_thresh(
            self,
            entrez_ids: List[int]
            ) -> Dict[str, int]:
        """Returns a dict of the phenotypes and number of hits that 
        are associated with the entrez ids. Uses thresh for each entrez id
        """
        dict_phenotype_to_hits = {}
        for phenotype in self.get_all_phenotypes():
            dict_phenotype_to_hits[phenotype] = 0
        for entrez_id in entrez_ids:
            for key, value in self.thresh_entrez_to_phenotypes(entrez_id).items():
                dict_phenotype_to_hits[key] += value
        return dict_phenotype_to_hits

if __name__ == "__main__":
    mpo = Gene_to_MPO_Features()
    while True:
        try:
            input_string = input("Enter entrez id: ")
        except (EOFError, KeyboardInterrupt) as e:
            print("Exiting...")
            exit()
        try:
            entrez_id = int(input_string)
            print("Phenotypes affected: "
                + str(mpo.entrez_to_phenotypes(entrez_id)))
        except ValueError:
            print("Invalid input")
        









