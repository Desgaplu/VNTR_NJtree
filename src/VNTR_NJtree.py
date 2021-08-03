#! python3
# -*- coding: utf-8
"""
Created on Mon Apr 27 22:39:15 2020

@author: pldesgagne

Description:
    Construct a Neighbor-Joining tree with VNTR data.

    Load VNTR data from a excel file.
    Scan and display information about the VNTR data.
    Calculate the genetic distance with either the Nei's distance formula
    or the Cavalli-Sforza chord distance formula.
    Saves a Neighbor-Joining tree that can be opened with a tree viewing
    program such as MEGA.

    The excel data should be formated as: (first line is ignored)
        Name    Locus   Allele # headers are ignored
        sample1 loci1   allele1
        sample1 loci2   allele1
        sample1 loci2   allele2
        sample2 loci1   allele1
"""
__version__ = "1.0.1"
import pandas as pd  # using the readExcel method
import tkinter as tk  # using the open filedialog
from tkinter import filedialog  # using the open filedialog
from math import pi  # for CavalliSforza algo
import numpy as np
import copy  # Copying np.matrix
from Bio.Phylo import write  # Biopython
from Bio.Phylo import BaseTree  # Biopython
import re

# Initialise tkinter to enable the uses of filedialog
root = tk.Tk()
root.withdraw()  # Prevent a empty window to be opened

# Title and icon for eventual GUI
# root.title("VNTR to Neighbor-Joining tree")
# root.iconbitmap('phylotree.ico')


class DistanceMatrix(object):
    """
    Distance matrix class.

    Contains a list of IDs and a distance matrix for the same IDs.
    IDs can be both indices or names (str).

    Attributes
    ----------
        names: a list of names (str)
    """

    def __init__(self, names):

        #  Initialize the matrix
        self.matrix = np.matrix(np.zeros(shape=(len(names),
                                                len(names)))).astype(float)
        self.names = names

        #  ID-to-index dictionary
        self.indices = {}  # {name:index}
        #  index-to-ID dictionary
        self.index_to_name = {}  # {index:name}

        # Create the indices dictionnary with the provided list of names
        self.__createIndices(names)

    def __createIndices(self, names):
        index = 0
        indices = {}
        index_to_name = {}
        for name in names:
            indices.update({name: index})
            index_to_name.update({index: name})
            index += 1
        self.indices = indices
        self.index_to_name = index_to_name

    #  Dict-like behaviour

    def __getitem__(self, item):
        """
        Get a distance between two sequences.

        Input
        -----
            item: a tuple of sequence names

        Return
        ------
            the genetic distance (float)
        """
        assert type(item) is tuple
        assert len(item) == 2
        # Verify if names are supplied instead of indices
        if all(isinstance(x, str) for x in item):
            return self.matrix[self.indices[item[0]], self.indices[item[1]]]
        else:
            return self.matrix[item[0], item[1]]

    def __setitem__(self, key, value):
        """
        Add an item to the matrix.

        Input
        -----
            a 2-item tuple
                key: name (str) or index (int)
                value: genetic distance (float)
        """
        assert type(key) is tuple
        assert len(key) == 2
        if all(isinstance(x, str) for x in key):
            self.matrix[self.indices[key[0]], self.indices[key[1]]] = value
            self.matrix[self.indices[key[1]], self.indices[key[0]]] = value
        else:
            self.matrix[key[0], key[1]] = value
            self.matrix[key[1], key[0]] = value

    def __delitem__(self, ids):
        """Remove a single ID."""
        self.matrix = np.delete(np.delete(self.matrix, ids, axis=0),
                                ids, axis=1)
        del self.names[ids]
        self.__createIndices(self.names)

    def __len__(self):
        """Return count of IDs in the matrix."""
        return len(self.indices)


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

        # DistanceMatrix
        self.distanceMatrix = None

        # Tree instance of Biopython Phylo BaseTree module
        self.tree = None

        # List of all loci names
        self.lociNames = ''
        self.lociCount = 0

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
        source_file_path = filedialog.askopenfilename(
            title="Select an Excel File",
            filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*"))
            )

        if source_file_path == '':
            raise CancelException("Open file cancelled")

        # Read the excel data into Pandas DataFrame
        fileData = pd.read_excel(source_file_path, header=0,
                                 names=["Name", "Locus", "Allele"],
                                 usecols="A:C",
                                 dtype={'Name': str,
                                        'Locus': str,
                                        'Allele': str})
        print("\nLoading: " + source_file_path)

        # load data into list of pop then into self.excelData
        pop_list = []
        # Initialize the progress bar at 0%
        self.__printProgressBar(0, fileData.shape[0], 'Loading Data:',
                                'Complete', 50)
        j = 0  # progress counter
        for i in fileData.itertuples(index=False):
            j += 1  # update progress counter
            for item in pop_list:
                if item.name == i.Name:
                    item.addLocus(i.Locus, i.Allele)
                    break
            else:
                currentPop = Pop(i.Name)
                currentPop.addLocus(i.Locus, i.Allele)
                pop_list.append(currentPop)
            # update the progress bar
            self.__printProgressBar(j, fileData.shape[0], 'Loading Data:',
                                    'Complete', 50)

        # save data into self.excelData
        self.excelData = pop_list
        print("\nExcel Data succesfully loaded:")

        # scan the data
        self.__scanExcelData(self.excelData)

    def __scanExcelData(self, data):
        """
        Print usefull information about the supplied VNTR data.

        Parse the excel data and show different information which may help the
        user find potential error in the provided excel files such as typos in
        names or wrong sample count.

        Input:
            data: the data to be scanned in format [Pop()]
        """
        lociNames = {}  # {key=Locus name: value=Locus count}
        numSamples = len(data)
        unusual_pop = []  # pop with less than minLociCount (potential typo)
        minLociCount = 3  # minimum loci to trigger a warning
        loci_missing_pop = []  # pop with missing loci
        non_int_allele = []  # "pop, locus, alelle" with non int allele value

        # Finds all the loci and their total count in the entire data
        for pop in data:
            pop_keys = pop.loci.keys()
            # saves loci with less than minLociCount (potential typo)
            if len(pop_keys) < minLociCount:
                unusual_pop.append(pop.name)
            for locus, alleles in pop.loci.items():
                if locus in lociNames:
                    lociNames[locus] += 1
                else:
                    lociNames[locus] = 1
                for allele in alleles.keys():
                    if not str(allele).isdecimal():
                        non_int_allele.append(
                            f'Sample {pop.name} at locus {locus} = "{allele}"')
        self.lociNames = lociNames.keys()
        self.lociCount = len(lociNames)

        # Print the number of samples found
        tabs = '\t'
        print(f"{tabs}{numSamples} samples were found.")

        # Print all the loci found and their count
        print(f"{tabs}{len(lociNames)} different loci were found:")
        tabs += '\t'
        print(f"{tabs}{'Locus name':<10}{'Count':>10}")
        for name, count in lociNames.items():
            if name is np.nan:
                print(f"{tabs}{'MISSING':<10}{count:>10}")
            else:
                print(f"{tabs}{name:<10}{count:>10}")

        # Print pop with one or more missing loci
        # Raise a error if a pop name or loci name is empty
        for pop in data:
            if pop.name is np.nan:
                print("\tERROR: One or more samples have no sample name.")
                raise CancelException("Add names or delete empty lines.")
            for locus in self.lociNames: 
                if locus is np.nan:
                    print("\tERROR: One or more samples have loci without a name.")
                    raise CancelException("Add names or delete empty lines.")
                if locus not in pop.loci:
                    loci_missing_pop.append(f"{pop.name} is missing locus " +
                                            locus)
        if loci_missing_pop:
            # If more that 80% of sample are missing a specific locus,
            # this locus may be a typo
            if len(loci_missing_pop) <= len(data)*0.8:
                print("\tWARNING:")
                for message in loci_missing_pop:
                    print(f"{tabs}{message}")
            else:
                print("\tWARNING: Check for locus name typos.")

        # Print pop with less than minLociCount (potential typo in pop name)
        if unusual_pop:
            print("\tWARNING: The following samples have less",
                  f"than {minLociCount} loci:")
            for pop in unusual_pop:
                print(f"{tabs}{pop}")

        # Print pop with non numerical allele values, raise error
        if non_int_allele:
            print("\tERROR: The following samples have alleles with",
                  "non integer values:")
            for allele in non_int_allele:
                print(f"{tabs}{allele}")
            raise CancelException("Fix the allele values.")

    def __printProgressBar(self, iteration, total, prefix='', suffix='',
                           length=100, fill='█', printEnd="\r"):
        """
        Call in a loop to create terminal progress bar.

        Input:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
        """
        percent = f"{100 * (iteration / total):.1f}"
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printEnd)
        # Print New Line on Complete
        if iteration == total:
            print()

    def saveTreeFile(self):
        """
        Save the tree data into a tree file.

        Ask the user for the name and destination of tree file created by
        this script.
        """
        print('\nEnter the destination for the tree file: ')
        dest_file_path = filedialog.asksaveasfilename(
            title="Select an destination and name for your Tree file",
            filetypes=(("Newick Tree", "*.nwk"), ("All files", "*.*"))
            )
        if dest_file_path == '':
            raise CancelException("Save file cancelled")
        # In case of file overwriting, prevent saving a ".nwk.nwk" file
        if dest_file_path[-4:] == '.nwk':
            dest_file_path = dest_file_path[:-4]
        write(njtree.tree, f"{dest_file_path}.nwk", 'newick')
        # Correct the apostrophes caharacter in final file
        # for MEGA software usage (otherwise it raise error)
        self.correct_newick(f"{dest_file_path}.nwk")
        print(f'\tTree saved in {dest_file_path}.nwk')

    def correct_newick(self, file_path):
        """
        Correct the Newick tree apostrophes.

        Correct the Newick tree apostrophes from "\'" to "''" to be readable
        with Mega (otherwise it raise a error). It overwrite the file with the
        modified version.

        Parameters
        ----------
        file_path : Str
            The file path of the Newick tree file (*.nwk)

        Returns
        -------
        None.

        """
        # open file
        with open(file_path, mode='r') as f:
            file = f.read()
            # replace the  apostrophes
            newfile = re.sub(r"\\'", "''", file)

        # overwrite the file
        with open(file_path, mode='w') as j:
            j.write(newfile)

    def buildTree(self, data=None, formula='Cavalli'):
        """
        Take the VNTR data and return a Neighbor-Joining tree.

        The distance matrix can be calculated with either the Cavalli-Sforza
        chord distance formula or the Nei's distance formula.

        Input:
            data (optionnal): the VNTR data to be analysed.
                Instance data is used by default.
            formula (optionnal): Cavalli (default) or Nei
        """
        # Verify if excel data is loaded / raise error otherwise
        if data is None and self.excelData:
            data = self.excelData
        else:
            raise AttributeError('No data was loaded into instance.')

        # Calculate the distance matrix from the loaded data according to
        # the specified formula
        algo = None
        if formula == 'Cavalli':
            algo = self.__cdCavalliSforza
        elif formula == 'Nei':
            algo = self.__neiDistance
        else:
            raise ValueError('"' + formula + '"' + " doesn't exist.")

        self.distanceMatrix = self.__geneticDistance(data, algo)

        # Build the tree from the distance matrix
        self.tree = self.__neighbor(self.distanceMatrix)

    def __cdCavalliSforza(self, dsum):
        """
        Return the Cavalli-Sforza chord distance.

        Uses the Cavalli-Sforza chord distance formula.
        Distance of 0 indicate that 2 sample are identical.
        Max distance is (2/pi)*sqrt(2) = 0.900316

        Input:
            dsum: sum of the squareroot of the multiplication of each allele
                frequency between 2 pop.

        Return
        ------
            Cavalli-Sforza chord distance
        """
        return (2/(pi*self.lociCount))*(2*abs(1-dsum))**0.5

    def __neiDistance(self, dsum):
        """
        Return Nei's DA distance.

        Uses the Nei's DA distance 1983 formula.
        Distance of 0 indicate that 2 sample are identical.
        Max distance is 1.

        Input:
            dsum: sum of the squareroot of the multiplication of each allele
                frequency between 2 pop.

        Return
        ------
            Nei's distance
        """
        return abs(1-dsum)/self.lociCount

    def __geneticDistance(self, data, algo):
        """
        Build a distance matrix with the VNTR data.

        Input:
            data : the VNTR data loaded by the loadExcelData method.
            algo : the formula used to calculate the genetic distance

        Returns
        -------
            DistanceMatrix

        """
        dmatrix = DistanceMatrix([pop.name for pop in data])  # Initialize
        dsum = 0  # sum of (allele frequency popA * allele frequency popB)**0.5
        distance = 0  # genetic distance between 2 pop
        num_pop = len(data)  # For visual progress
        current_pop = 1  # For progress bar

        # Initialize the progress bar at 0%
        print()
        self.__printProgressBar(0, num_pop, 'Calculating Distances:',
                                'Complete', 50)
        # for each pop
        for pop in data:
            self.__printProgressBar(current_pop, num_pop,
                                    'Calculating Distances:',
                                    f'Complete ({current_pop}/{num_pop})',
                                    50)
            current_pop += 1
            # compare to all previous samples excluding self
            # for versus in data:                   # Square matrix
            for versus in data[:data.index(pop)]:   # Lower triangular matrix
                distance = 0  # initialize distance for pop vs versus
                remainingLocus = self.lociCount
                # for each locus of pop
                for locus in pop.loci:
                    remainingLocus -= 1
                    dsum = 0  # initialize sum for a single locus
                    # for each allele in locus of pop
                    for allele in pop.loci[locus]:
                        # Check if versus has current pop locus
                        if locus in versus.loci:
                            # Check if versus has same allele
                            if allele in versus.loci[locus]:
                                # sum of sqrt of allele frequencies product
                                dsum += (pop.frequency(locus)[allele] *
                                         versus.frequency(locus)[allele])**0.5
                    # Adding the distance for this locus
                    distance += algo(dsum)  # use the supplied formula fonction
                # Adding distance for missing locus
                distance += algo(0) * remainingLocus
                # Place calculated distance between pop and versus in matrix
                dmatrix[pop.name, versus.name] = distance
        return dmatrix

    def __neighbor(self, distance_matrix):
        """
        Construct and return a Neighbor-Joining tree.

        Input:
            distance_matrix : a DistanceMatrix instance

        Returns
        -------
            Bio.Phylo.BaseTree instance
        """
        print("\nStarting Neighbor-Joining.")
        # Formulas for the neighbor-joining matrix and minimum pair
        rptsum = lambda arr: np.repeat(np.sum(arr)/(np.size(arr)-2),
                                       np.size(arr))
        mapvsum = lambda mat: np.matrix([rptsum(line) for line in mat])
        idxmin = lambda mat: np.unravel_index(np.argmin(mat), np.shape(mat))

        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        tot_len = len(distance_matrix)  # for progress bar

        # init terminal clades
        clades = [BaseTree.Clade(None, name) for name in dm.names]

        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        # special cases for Minimum Alignment Matrices
        if len(dm) == 1:
            root = clades[0]

            return BaseTree.Tree(root, rooted=False)
        elif len(dm) == 2:
            # minimum distance will always be [1,0]
            min_i = 1
            min_j = 0
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            clade1.branch_length = dm[min_i, min_j] / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length
            inner_clade = BaseTree.Clade(None, "Inner")
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            clades[0] = inner_clade
            root = clades[0]

            return BaseTree.Tree(root, rooted=False)

        # Initialize the progress bar at 0%
        self.__printProgressBar(0, tot_len, 'Joining:', 'Complete', 50)
        while len(dm) > 2:
            # Progress bar update
            current_pos = tot_len - len(dm)
            self.__printProgressBar(current_pos+3, tot_len, 'Joining:',
                                    f'Complete ({current_pos+3}/{tot_len})',
                                    50)

            # calculate prerequisites for neighbor-joining matrix
            SH = mapvsum(dm.matrix)
            SV = SH.transpose()

            # Build the neighbor-joining matrix M
            Id = np.identity(len(dm.matrix))
            M = dm.matrix + (np.multiply(Id, SH + SV) - SH - SV)

            # Find minimum distance pair
            min_i, min_j = idxmin(M+Id)  # +Id to prevent min_i == min_j

            # create clades with the minimum distance pair found
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = BaseTree.Clade(None, "Inner" + str(inner_count))
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)

            # assign branch lengths
            clade1.branch_length = (
                dm[min_i, min_j] + SH[min_i, min_j] - SV[min_i, min_j]
            ) / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length

            # update clades list with new clade pair
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new clade at the index of min_j
            u = [(dm[min_i, k] + dm[min_j, k] - dm[min_i, min_j]) / 2 for k in
                 range(len(dm))]
            dm.matrix[min_j] = u
            dm.matrix[:, min_j] = np.matrix(u).transpose()
            dm.names[min_j] = "Inner" + str(inner_count)
            del dm[min_i]

        # set the last clade as one of the child of the inner_clade
        root = None
        if clades[0] == inner_clade:
            clades[0].branch_length = 0
            clades[1].branch_length = dm[1, 0]
            clades[0].clades.append(clades[1])
            root = clades[0]
        else:
            clades[0].branch_length = dm[1, 0]
            clades[1].branch_length = 0
            clades[1].clades.append(clades[0])
            root = clades[1]

        return BaseTree.Tree(root, rooted=False)

    def alleleFrequency(self):
        """Print an allele frequency table."""
        for pop in self.excelData:
            for locus in pop.loci:
                for allele, freq in pop.frequency(locus).items():
                    print(pop.name, locus, allele, round(freq, ndigits=4))

    def executeCommand(self):
        """Excute the workflow."""
        try:
            # Load an excel file contain VNTR data
            self.loadExcelData()
            query = input('\nIs the displayed information correct? [y/n] ')
            if query.lower() != 'y':
                raise CancelException(
                    "VNTR information considered incorrect by user.")
            else:
                # Build a phylogenetic tree
                # self.buildTree(formula='Nei')
                self.buildTree(formula='Cavalli')
                print('\nNeighbor-Joining tree constructed.')
                # Save the tree in specified file
                self.saveTreeFile()
                input('\nPress Enter to exit.')
        except CancelException as e:
            print(f'\n***{e.message}***')
            input('Press Enter to exit.')
        except Exception as e:
            print(f'\n***Unexpected error: {e}***')
            input('Press Enter to exit.')


class Pop():
    """
    Representation of a single sample.

    Attributes
    ----------
            name: Name of sample (str)
    """

    def __init__(self, name):

        self.name = name

        # loci: a dict of each locus name containing alleles values.
        # {str_loci_name : {int_allele_value : int_count, }, }
        self.loci = {}

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


class CancelException(BaseException):
    """Exception for dealing with the cancel command from user."""

    def __init__(self, message):
        self.message = message


if __name__ == '__main__':

    njtree = NJTreeConstructor()
    njtree.executeCommand()
