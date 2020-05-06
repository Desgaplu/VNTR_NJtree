#! python3
"""
Created on Mon Apr 27 22:39:15 2020

@author: pldesgagne

Description:
    Blah # TODO
"""
# TODO add pandas and tkinter in UML diag
import pandas as pd
import tkinter as tk
from tkinter import filedialog
# from Bio import Phylo as phy

# Initialise tkinter to enable the uses of filedialog
root = tk.Tk()
root.withdraw()  # Prevent a empty window to be opened

# Title and icon for eventual GUI
# root.title("VNTR to Neighbor-Joining tree")
# root.iconbitmap('phylotree.ico')


class NJTreeConstructor():
    """
    Construct a Neighbor-Joining tree with VNTR data.

    Creates a object that can load VNTR data from a excel file and return
    a Neighbor-Joining tree that can be opened with a tree viewing program
    such as MEGA.

    The excel data should be formated as: (first line is ignored)
        Name    Locus   Allele # headers are ignored
        sample1 loci1   allele1
        sample1 loci2   allele1
        sample1 loci2   allele2
        sample2 loci1   allele1
        ...
    """

    def __init__(self):

        # VNTR data in a list of Pop()
        self.excelData = None

        # DistanceMatrix from the Biopython Phylo TreeConstruction module
        self.distanceMatrix = None

        # Tree from the Biopython Phylo BaseTree module
        self.tree = None

    def loadExcelData(self):
        """
        Load the data from excel sheet into a list of Pop().

        The excel sheet contains name, locus and alleles for each samples.
        Saves results into class variable self.excelData

        Input:
            source: file name and address for the file containing the samples
                name and their allele value for each locus.
                format in each colomn should be:
                    Name Locus Allele # headers names are ignored
                    name1 loci1 allele1
                    name1 loci2 allele2
        """
        print("Select the excel file containing the VNTR data: ")
        print('(The excel file should contain',
              '3 columns "Name", "Locus", "Allele".)')

        # Ask for the excel file path
# =============================================================================
#         source_file_path = filedialog.askopenfilename(
#             title="Select an Excel File",
#             filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*"))
#             )
# =============================================================================
        # TEMPORARY FOR TESTING, use above ^ # TODO
        source_file_path = 'data.xlsx'

        # Read the excel data into Pandas DataFrame
        fileData = pd.read_excel(source_file_path, header=0,
                                 names=["Name", "Locus", "Allele"],
                                 usecols="A:C",
                                 dtype={'Name': str,
                                        'Locus': str,
                                        'Allele': int})
        print("\nLoading: " + source_file_path)

        # load data into list of pop then into self.excelData
        pop_list = []
        for i in fileData.itertuples(index=False):
            for item in pop_list:
                if item.name == i.Name:
                    item.addLocus(i.Locus, i.Allele)
                    break
            else:
                currentPop = Pop(i.Name)
                currentPop.addLocus(i.Locus, i.Allele)
                pop_list.append(currentPop)

        # save data self.excelData
        self.excelData = pop_list
        print("Excel Data succesfully loaded:")

        # scan the data
        self.__scanExcelData(self.excelData)

    def __scanExcelData(self, data):
        """
        Print usefull information about the supplied VNTR data.

        Parse the excel data and show different information which may help the
        user find potential error in the provided excel files such as typos in
        names or wrong sample number.

        Input:
            data: the data to be scanned in format [Pop()]
        """
        # print a data feedback (number of samples, loci names etc)
        lociNames = {}  # {key=Locus name: value=Locus count}
        numSamples = len(data)
        unusual_pop = []  # pop with less than minLociCount (potential typo)
        minLociCount = 6

        # Finds all the Loci and their count
        for pop in data:
            pop_keys = pop.loci.keys()
            if len(pop_keys) < minLociCount:
                unusual_pop.append(pop.name)
            for locus in pop_keys:
                if locus in lociNames:
                    lociNames[locus] += 1
                else:
                    lociNames[locus] = 1

        # Print the number of samples found
        tabs = '\t'
        print(f"{tabs}{numSamples} samples were found.")

        # Print all the loci found and their count
        print(f"{tabs}{len(lociNames)} different loci were found:")
        tabs += '\t'
        print(f"{tabs}Locus name\t\t Count")
        for name, count in lociNames.items():
            print(f"{tabs}{name:<8}{count:>10}")

        # Print pop with less than minLociCount (potential typo in pop name)
        if unusual_pop:
            print(f"\tWARNING: The following samples have less",
                  "than {minLociCount} loci:")
            for locus in unusual_pop:
                print(f"{tabs}{locus}")

    def saveTreeFile(self, destination: str):
        """
        Save the tree data into a tree file.

        Input:
            destination: file name and address for the tree file created by
                this script.
        """
        pass  # TODO

    def buildTree(self, data=None):
        """
        Take the VNTR data and return a Neighbor-Joining tree.

        Uses the Cavalli-Sforza chord distance formula to create a distance
        matrix then the Biopython Phylo package to build the neighbor tree

        Input
        -----
            data (optionnal): the VNTR data to be analysed. If none is passed,
                the self.excelData from the instance is used.

        Return
        ------
            Tree from the Biopython Phylo BaseTree module.
        """
        # verify if excel data is loaded / raise error otherwise
        if data is None and self.excelData:
            data = self.excelData
        else:
            raise AttributeError('No data was loaded into instance.')

        # Calculate the distance matrix from the loaded data
        self.distanceMatrix = self.__cdCavalliSforza(data)

        # Build the tree from the distance matrix
        self.tree = self.__neighbor(self.distanceMatrix)

    def __cdCavalliSforza(self, data):
        """
        Build a distance matrix with the VNTR data.

        Uses the Cavalli-Sforza chord distance formula.

        Input
        -----
            data: the VNTR data loaded by the loadExcelData method.

        Return
        ------
            DistanceMatrix from the Biopython Phylo TreeConstruction module
        """
        pass  # TODO

    def __neighbor(self, matrix):
        """
        Apply the neighbor joining algo on a distance matrix.

        Input
        -----
            matrix : DistanceMatrix from the Biopython Phylo TreeConstruction
                module

        Returns
        -------
            Tree from the Biopython Phylo BaseTree module
        """
        pass  # TODO


class Pop():
    """
    Representation of a single sample.

    Attributes
    ----------
            name: Name of sample in string format
            loci: a dict of each locus name containing alleles values.
                {loci_name:{allele_value:count,},}
    """

    def __init__(self, name):

        self.name = name
        self.loci = {}  # dict {name:{allele:count,},}

    def addLocus(self, locus: str, allele: int):
        """
        Add a locus and its allele value in the loci dict.

        If the locus is already present, the new allele is added to the dict
        If the allele is already present for that loci, +1 to count.
        """
        if locus not in self.loci:
            self.loci[locus] = {allele: 1}
        elif allele in self.loci[locus]:
            self.loci[locus][allele] += 1
        else:
            self.loci[locus][allele] = 1

    def frequency(self, locus):
        """
        Return the allelic frequency.

        Return a dict {allele:frequency} of the frequency of each allele
        for the specified locus.

        Ex.: for 3 alleles {650:2, 720:1} => {650:0.6666, 720:0.3333}
        (Divides the count of each allele by the total number of allele for
        the specified loci.)
        """
        nb_allele = 0
        for allele in self.loci[locus]:
            nb_allele += self.loci[locus][allele]
        return {i: (self.loci[locus][i]/nb_allele) for i in self.loci[locus]}

    def __str__(self):
        """Return the name of the pop."""
        return self.name


def testUnit():
    """
    For testing the script.

    Returns
    -------
    None.

    """
    sample = Pop('R1105')
    sample.addLocus('Tr1', 650)
    sample.addLocus('Tr1', 650)
    sample.addLocus('Tr1', 720)
    sample.addLocus('Tr2', 150)
    # assert test_bool, "Test failed comment"
    # print("All test passed")


if __name__ == '__main__':

    njtree = NJTreeConstructor()
    njtree.loadExcelData()
# =============================================================================
#     query = input('Is the displayed information correct? [y/n] ')
#     if query.lower() != 'y':
#         print('\nCancelling...')
#         input('Press Enter to exit.')
#     else:
#         njtree.buildTree()
#         print('\nNeighbor-Joining tree constructed.')
#         destination = input('Enter a name for the tree file: ')
#         njtree.saveTreeFile(destination)
#         print(f'\nTree saved in {destination}')
#         input('Press Enter to exit.')
# =============================================================================
