"""Gene to Pathways V0.2

=== TCAG Summer Project 2018 ===
Yoonsik Park
BCB330Y1
University of Toronto

=== Module Description ===
Given a certain human entrez ID, the module will return all pathways.
This module also supports returning pathways for multiple entrez ids,
and also supports filtering the pathways based on gene numbers
"""

from typing import Dict, Tuple, List, Optional
import csv
import linecache

PATHWAY_FILENAME = "allSizes_GOincludingIEA_pathways_20180213.GMT"
DELIMITER = "\t"

class Gene_to_Pathways:
    """This class loads the pathways file efficiently and allows for quick
    lookups of entrez genes to pathways

    === Attributes ===
    filename:
        The filename containing the pathway data. Most likely:
        "allSizes_GOincludingIEA_pathways_20180213.GMT"
    delimiter:
        The delimiter to be used. (i.e. comma ',' or tab '\t')

    === Private attributes ===
    _dict_genes_to_lines:
        A dictionary to map entrez gene IDs to the corresponding line numbers 
        of the file.
        i.e. "23" --> [5, 6, 19]
    _list_number_genes_by_line:
        A list storing the number of genes for each pathway, ordered by line 
        number.
        _list_number_genes_by_line[i] = number of genes on line number (i + 1)
    _dict_databases_to_lines:
        A dictionary to map databases supported by the file to a set of line 
        numbers
        i.e. "GO", "KEGG", "NCI" --> set([1, 3, 4])
    """
    filename: str
    delimiter: str
    _dict_genes_to_lines: Dict[int, List[int]]
    _list_number_genes_by_line: List[int]
    _dict_databases_to_lines: Dict[str, set]

    def __init__(
            self,
            filename: str = PATHWAY_FILENAME,
            delimiter: str = DELIMITER
            ) -> None:
        """Opens the pathway file, and creates data structures based on it."""
        self.filename = filename
        self.delimiter = delimiter
        self._dict_genes_to_lines = {}
        self._list_number_genes_by_line = []
        self._dict_databases_to_lines = {}
        file = open(filename, "r")
        file_reader = csv.reader(file, delimiter = self.delimiter)
        # firstline = True
        line_number = 0
        for row in file_reader:
            line_number += 1
            # if firstline:
            #     firstline = False
            #     continue
            gene_set = row[0].split(':')[0].upper()
            if gene_set in self._dict_databases_to_lines:
                self._dict_databases_to_lines[gene_set].add(line_number)
            else:
                self._dict_databases_to_lines[gene_set] = set([line_number])
            
            self._list_number_genes_by_line.append(len(row) - 2)

            for str_entrez_id in row[2:]:
                entrez_id = int(str_entrez_id)
                if entrez_id in self._dict_genes_to_lines:
                    self._dict_genes_to_lines[entrez_id].append(line_number)
                else:
                    self._dict_genes_to_lines[entrez_id] = [line_number]
        file.close()

    def get_all_databases(self) -> List[str]:
        """Returns the list of databases found in the file."""
        return sorted(self._dict_databases_to_lines)

    def line_num_description(
            self,
            line_num: int,
            delimiter: str = ' '
            ) -> str:
        """Returns a description of the gene set on the corresponding
        line number
        """
        row = linecache.getline(self.filename, line_num).split(self.delimiter)
        return row[0] + delimiter + row[1]

    def _entrez_to_line_nums(
            self,
            entrez_id: int,
            database: str,
            filter_min: int = 1,
            filter_max: int = 50000
            ) -> List[int]:
        """ Returns a list of line numbers corresponding to any gene sets that
        are associated with the entrez id. The database parameter is required
        to filter based on gene set database. Optionally, filter_min and 
        filter_max are used to filter the gene sets such that the number of 
        genes inside a gene set are within the range (inclusive). Returns 
        an empty list if the entrez id doesn't exist or no gene sets are found 
        that match the criteria.
        """
        database = database.strip().upper()
        line_nums = []

        if entrez_id not in self._dict_genes_to_lines:
            return line_nums
        if database not in self._dict_databases_to_lines:
            print(__name__ + ": Warning '" + database + "' database not found!")
            return line_nums

        for line_num in self._dict_genes_to_lines[entrez_id]:
            num_genes = self._list_number_genes_by_line[line_num - 1]
            if filter_min <= num_genes <= filter_max and \
                    line_num in self._dict_databases_to_lines[database]:
                line_nums.append(line_num)
        return line_nums

    def entrez_to_num_gene_sets(
            self,
            entrez_id: int,
            database: str,
            filter_min: int = 1,
            filter_max: int = 50000
            ) -> int:
        """ Returns the number of gene sets that are associated with the 
        entrez id. The database parameter is required to filter based on gene 
        set database. Optionally, filter_min and filter_max are used to 
        filter the gene sets such that the number of genes inside a gene 
        set are within the range (inclusive). Returns 0 if the entrez id 
        doesn't exist or no gene sets are found that match the criteria.
        """
        return len(self._entrez_to_line_nums(entrez_id, database,
                                            filter_min = filter_min,
                                            filter_max = filter_max))

    def entrez_to_gene_sets(
            self,
            entrez_id: int,
            database: str,
            filter_min: int = 1,
            filter_max: int = 50000
            ) -> List[Tuple[int, str]]:
        """Returns a list of tuples of the gene sets that are associated with 
        the entrez id. Each tuple contains the line number and gene set 
        description respectively. The database parameter is required to filter 
        based on gene set database. Optionally, filter_min and filter_max are 
        used to filter the gene sets such that the number of genes inside a 
        gene set are within the range (inclusive). Returns an empty list if 
        the entrez id doesn't exist or no gene sets are found that match the 
        criteria.
        """
        
        gene_sets = []
        for line_num in self._entrez_to_line_nums(entrez_id, database,
                                                filter_min = filter_min,
                                                filter_max = filter_max):
            description = self.line_num_description(line_num, delimiter = ';')
            gene_sets.append((line_num, description))
        return gene_sets

    def multi_entrez_to_gene_sets(
            self,
            entrez_ids: List[int],
            database: str,
            filter_min: int = 1,
            filter_max: int = 30000
            ) -> List[Tuple[int, str, int]]:
        """Returns a list of tuples of the gene sets and number of hits that 
        are associated with the entrez ids. Each tuple contains the line 
        number, gene set description, and counts/hits respectively. The 
        database parameter is required to filter based on gene set database. 
        Optionally, filter_min and filter_max are used to filter the gene sets 
        such that the number of genes inside a gene set are within the 
        range (inclusive). Returns an empty list if no gene sets are found or 
        no gene sets are found that match the criteria.
        """

        # dict_lines_to_hits maps a line number (i.e. the gene set)
        # to the number of hits it has received
        dict_lines_to_hits = {}

        for entrez_id in entrez_ids:
            for line_num in self._entrez_to_line_nums(entrez_id, database,
                                                    filter_min = filter_min,
                                                    filter_max = filter_max):
                if line_num in dict_lines_to_hits:
                    dict_lines_to_hits[line_num] += 1
                else:
                    dict_lines_to_hits[line_num] = 1

        list_return = []
        for line_num, hits in dict_lines_to_hits.items():
            description = self.line_num_description(line_num, delimiter = ';')
            list_return.append((line_num, description, hits))

        return list_return

if __name__ == "__main__":
    pathways = Gene_to_Pathways()
    print("Databases found: " + str(pathways.get_all_databases()))
    while True:
        try:
            db_string = input("Enter database: ")
            input_string = input("Enter entrez id: ")
        except (EOFError, KeyboardInterrupt) as e:
            print("Exiting...")
            exit()
        try:
            entrez_id = int(input_string)
            print("Pathways affected: "
                  + str(pathways.entrez_to_gene_sets(entrez_id, db_string)))
        except ValueError:
            print("Invalid input")
        









