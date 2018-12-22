"""Gene Interval to Repetitive Element Features V0.4

=== TCAG Summer Project 2018 ===
Yoonsik Park
BCB330Y1
University of Toronto

=== Module Description ===
Given a certain chromosome, start, and end location, this module will return
repetitive element features. The repetitive element file is by default 
assumed to be compressed.
"""

from typing import Dict, Tuple, List, Optional, Set
from intervaltree import Interval, IntervalTree
import csv
import linecache
import gzip

REPEAT_ELEM_FILENAME = "RLCRs_DNN-CNV.txt.gz"
DELIMITER = "\t"
# Column numbers starting from zero
CHROM_COLUMN = 0
START_COLUMN = 1
END_COLUMN = 2
REPEAT_ELEM_COLUMN = 3

class Gene_Interval_to_Repetitive_Elements:
    """This class loads the repetitive elements file and allows for
    quick lookups of matching intervals

    === Attributes ===
    filename:
        The filename containing the Repetitive Element data. Most likely: 
        "RLCRs_DNN-CNV.txt.gz"
    delimiter:
        The delimiter to be used. (i.e. comma or tab)

    === Private attributes ===
    _dict_chrom_to_trees:
        A dictionary where the chromosome (e.g. 'chrX')
        maps to an interval tree. Each interval in the tree contains the
        the repetitive element string
    """
    filename: str
    cache_filename: str
    delimiter: str
    _dict_chrom_to_trees: Dict[str, IntervalTree]
    set_of_element_types: Set[str]

    def __init__(
            self,
            filename: str = REPEAT_ELEM_FILENAME,
            delimiter: str = DELIMITER
            ) -> None:
        """Opens the Repetitive elements file, and creates data structures based on it."""
        self.filename = filename
        self.delimiter = delimiter
        self._dict_chrom_to_lines = {}
        self.set_of_element_types = set()
        # if the file is not compressed change to this
        # file = open(filename, "r")
        file = gzip.open(filename, "rt")
        file_reader = csv.reader(file, delimiter = delimiter)
        firstline = True
        line_number = 0
        for row in file_reader:
            line_number += 1
            if line_number % 1000000 == 0:
                print(__name__ + ": Importing line " + str(line_number))
            if firstline:
                firstline = False
                continue
            chrom = row[CHROM_COLUMN]
            start = int(row[START_COLUMN])
            end = int(row[END_COLUMN])
            description = row[REPEAT_ELEM_COLUMN]
            self.set_of_element_types.add(description)
            if chrom not in self._dict_chrom_to_lines:
                self._dict_chrom_to_lines[chrom] = IntervalTree()
            self._dict_chrom_to_lines[chrom][start:end] = description
        file.close()
        

    def get_all_intersecting_elements(self, chrom: str, start: int, end: int) -> Dict[str, int]:
        """Returns a dictionary of all the repetitive elements that intersect on chromosome
        <chrom>, and corresponding <start> and <end>. The dictionary maps the
        repetitive element name (e.g. "LINE") to the count (e.g. 3). Returns an empty list
        if invalid
        """
        dict_elements = {}
        for element in self.set_of_element_types:
            dict_elements[element] = 0
        if start > end or chrom not in self._dict_chrom_to_lines:
            return dict_elements

        for interval in self._dict_chrom_to_lines[chrom].search(start, end):
            dict_elements[interval.data] += 1
        return dict_elements
        







