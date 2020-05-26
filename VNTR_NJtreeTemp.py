#! python3
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
import tkinter as tk
from tkinter import filedialog
from math import pi
import numpy as np
import copy
from Bio.Phylo import write
from Bio.Phylo import BaseTree

# Initialise tkinter to enable the uses of filedialog
root = tk.Tk()
root.withdraw()  # Prevent a empty window to be opened

# Title and icon for eventual GUI
# root.title("VNTR to Neighbor-Joining tree")
# root.iconbitmap('phylotree.ico')

class DistanceMatrix(object):
    '''
    Distance matrix class.
    Contains a list of IDs and a distance matrix for the same IDs. Can produce
    a submatrix when given an ID list, or calculate ID lists of a desired size
    using DJ algorithm
    '''
    def __init__(self, names):
        #  Initialize the matrix
        self.matrix = np.matrix(np.zeros(shape=(len(names), len(names)))).astype(float)
        self.names = names
        
        #  ID-to-index dictionary
        self.indices = {} #  {name:index}
        #  index-to-ID dictionary
        self.index_to_name = {} #  {index:name}
        
        self.__createIndices(names)
        
    def __createIndices(self, names):
        index = 0
        indices = {}
        index_to_name = {}
        for name in names:
            indices.update({name:index})
            index_to_name.update({index:name})
            index += 1
        self.indices = indices
        self.index_to_name = index_to_name
    
    #  Dict-like behaviour

    def __getitem__(self, item):
        '''
        Get a distance between two sequences
        :param item: a tuple of sequence names
        :return:
        '''
        assert type(item) is tuple
        assert len(item) == 2
        # Verify if names are supplied instead of indices
        if all(isinstance(x, str) for x in item):
            return self.matrix[self.indices[item[0]], self.indices[item[1]]]
        else:
            return self.matrix[item[0], item[1]]

    def __setitem__(self, key, value):
        '''
        Add an item to the matrix
        :param key: a 2-item tuple
        :param item: float
        '''
        assert type(key) is tuple
        assert len(key) == 2
        if all(isinstance(x, str) for x in key):
            self.matrix[self.indices[key[0]], self.indices[key[1]]] = value
            self.matrix[self.indices[key[1]], self.indices[key[0]]] = value
        else:
            self.matrix[key[0], key[1]] = value
            self.matrix[key[1], key[0]] = value
            
    def __delitem__(self, ids):
        self.matrix = np.delete(np.delete(self.matrix, ids, axis=0), ids, axis=1)
        del self.names[ids]
        self.__createIndices(self.names)
        
    def __len__(self):
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

        # Tree from the Biopython Phylo BaseTree module
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
        for i in fileData.itertuples(index=False):
            for item in pop_list:
                if item.name == i.Name:
                    item.addLocus(i.Locus, i.Allele)
                    break
            else:
                currentPop = Pop(i.Name)
                currentPop.addLocus(i.Locus, i.Allele)
                pop_list.append(currentPop)
            print('Loading pop')

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
        names or wrong sample number.

        Input:
            data: the data to be scanned in format [Pop()]
        """
        lociNames = {}  # {key=Locus name: value=Locus count}
        numSamples = len(data)
        unusual_pop = []  # pop with less than minLociCount (potential typo)
        minLociCount = 4  # minimum loci to trigger a warning
        loci_missing_pop = []  # pop with missing loci
        non_int_allele = [] # "pop, locus, alelle" with non int allele value

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
                    if not allele.isdecimal():
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
        print(f"{tabs}Locus name\t\t Count")
        for name, count in lociNames.items():
            print(f"{tabs}{name:<8}{count:>10}")

        # Print pop with one or more missing loci
        for pop in data:
            for locus in self.lociNames:
                if locus not in pop.loci:
                    loci_missing_pop.append(f"{pop.name} is missing locus " +
                                            locus)
        if loci_missing_pop:
            # If more that 80% are missing a locus, this locus may be a typo
            if len(loci_missing_pop) <= len(data)*0.8:
                print("\tWARNING:")
                for message in loci_missing_pop:
                    print(f"{tabs}{message}")
            else:
                print("\tWARNING: Check for locus name typos.")

        # Print pop with less than minLociCount (potential typo in pop name)
        if unusual_pop:
            print(f"\tWARNING: The following samples have less",
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

    def saveTreeFile(self):
        """
        Save the tree data into a tree file.

        Input:
            destination: file name and address for the tree file created by
                this script.
        """
        print('\nEnter the destination for the tree file: ')
        dest_file_path = filedialog.asksaveasfilename(
            title="Select an destination and name for your Tree file",
            filetypes=(("Newick Tree", "*.nwk"), ("All files", "*.*"))
            )
        if dest_file_path == '':
            raise CancelException("Save file cancelled")
        if dest_file_path[-4:] == '.nwk':  # prevent saving a ".nwk.nwk" file
            dest_file_path = dest_file_path[:-4]
        write(njtree.tree, f"{dest_file_path}.nwk", 'newick')
        print(f'\tTree saved in {dest_file_path}.nwk')

    def buildTree(self, data=None, formula='Cavalli'):
        """
        Take the VNTR data and return a Neighbor-Joining tree.

        Uses the Cavalli-Sforza chord distance formula to create a distance
        matrix then the Biopython Phylo package to build the neighbor tree

        Input:
            data (optionnal): the VNTR data to be analysed. If none is passed,
                the self.excelData from the instance is used.

        Return
        ------
            Tree from the Biopython Phylo BaseTree module.
        """
        # Verify if excel data is loaded / raise error otherwise
        if data is None and self.excelData:
            data = self.excelData
        else:
            raise AttributeError('No data was loaded into instance.')

        # Calculate the distance matrix from the loaded data according to
        # the asked formula
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
        DEBUG = False  # Print an output for debugging
        num_pop = len(data) # For visual progress
        current_pop = 1 # For visual progress

        def printx(*x):
            if DEBUG:
                print(*x)

        # for each pop
        for pop in data:
            printx('Current: ', pop)
            print(f"Analysing Pop #{current_pop} of {num_pop}.({int((current_pop/num_pop)*100)}%)")
            current_pop += 1
            # compare to all previous samples excluding self
            # for versus in data:                   # Square matrix
            for versus in data[:data.index(pop)]:   # Lower triangular matrix
                distance = 0  # initialize distance for pop vs versus
                remainingLocus = self.lociCount
                printx('\tversus:', versus)
                # for each locus of pop
                for locus in pop.loci:
                    remainingLocus -= 1
                    dsum = 0  # initialize sum for a single locus
                    printx('\t'*2, locus)
                    # for each allele in locus of pop
                    for allele in pop.loci[locus]:
                        printx('\t'*3, allele)
                        # Check if versus has current pop locus
                        if locus in versus.loci:
                            # Check if versus has same allele
                            if allele in versus.loci[locus]:
                                printx('\t'*4, "present in", versus)
                                # sum of sqrt of allele frequencies product
                                dsum += (pop.frequency(locus)[allele] *
                                         versus.frequency(locus)[allele])**0.5
                            else:
                                printx('\t'*4, "absent in", versus)
                        else:
                            printx('\t'*4, versus, "doesnt have a", locus)
                    # Adding the distance for this locus
                    distance += algo(dsum)  # use the supplied formula fonction
                # Adding distance for missing locus
                distance += algo(0) * remainingLocus
                # Place calculated distance between pop and versus in matrix
                dmatrix[pop.name, versus.name] = distance
            printx('---------')
        printx(dmatrix)
        printx('---------')
        return dmatrix

    def __neighbor(self, matrix):
        """
        Apply the neighbor joining algo on a distance matrix.

        Input
        -----
            matrix : DistanceMatrix
                module

        Returns
        -------
            Tree instance of Bio.Phylo.BaseTree
        """
        
        print("Starting Neignbor.")
        tree = self.__nj(matrix)
        print("Neighbor ended.")
        return tree

    def __nj(self, distance_matrix): # Second attempt at fast and correct
        """Construct and return a Neighbor Joining tree.

        Input
        -----
            distance_matrix : a distance matrix

        Returns
        -------
            Tree from the Biopython Phylo BaseTree module
        """
        rptsum  = lambda arr: np.repeat(np.sum(arr)/(np.size(arr)-2),
                                        np.size(arr))
        mapvsum = lambda mat: np.matrix([rptsum(line) for line in mat])
        idxmin  = lambda mat: np.unravel_index(np.argmin(mat), np.shape(mat))
        
        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        
        # init terminal clades
        clades = [BaseTree.Clade(None, name) for name in dm.names]
        
        # init node distance
        #node_dist = [0] * len(dm)
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
        
        
        while len(dm) > 2:
            print("Joining pop #", len(dm))
            # calculate nodeDist
            # for i in range(0, len(dm)):
            #     node_dist[i] = 0
            #     for j in range(0, len(dm)):
            #         node_dist[i] += dm[i, j]
            #     node_dist[i] = node_dist[i] / (len(dm) - 2)
            SH = mapvsum(dm.matrix)
            SV = SH.transpose()
            
            # find minimum distance pair
            # min_dist = dm[1, 0] - node_dist[1] - node_dist[0]
            # min_i = 0
            # min_j = 1
            # for i in range(1, len(dm)):
            #     for j in range(0, i):
            #         temp = dm[i, j] - node_dist[i] - node_dist[j]
            #         if min_dist > temp:
            #             min_dist = temp
            #             min_i = i
            #             min_j = j
            I = np.identity(len(dm.matrix))
            M = dm.matrix + (np.multiply(I, SH + SV) - SH - SV)
            min_i, min_j = idxmin(M)
            # create clade
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = BaseTree.Clade(None, "Inner" + str(inner_count))
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            # assign branch length
            clade1.branch_length = (
                dm[min_i, min_j] + SH[min_i,min_j] - SV[min_i,min_j]
            ) / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            # for k in range(0, len(dm)):
            #     if k != min_i and k != min_j:
            #         dm[min_j, k] = (
            #             dm[min_i, k] + dm[min_j, k] - dm[min_i, min_j]
            #         ) / 2.0
            u = [(dm[min_i,k] + dm[min_j,k] - dm[min_i,min_j]) / 2 for k in 
                 range(len(dm))]
            dm.matrix[min_j] = u
            dm.matrix[:,min_j] = np.matrix(u).transpose()
            
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
    
# =============================================================================
#     def __nj2(self, distance_matrix): # Functionnal and correct but slow
#         """Construct and return a Neighbor Joining tree.
# 
#         :Parameters:
#             distance_matrix : DistanceMatrix
#                 The distance matrix for tree construction.
# 
#         """
#         # make a copy of the distance matrix to be used
#         dm = copy.deepcopy(distance_matrix)
#         # init terminal clades
#         clades = [BaseTree.Clade(None, name) for name in dm.names]
#         # init node distance
#         node_dist = [0] * len(dm)
#         # init minimum index
#         min_i = 0
#         min_j = 0
#         inner_count = 0
#         # special cases for Minimum Alignment Matrices
#         if len(dm) == 1:
#             root = clades[0]
# 
#             return BaseTree.Tree(root, rooted=False)
#         elif len(dm) == 2:
#             # minimum distance will always be [1,0]
#             min_i = 1
#             min_j = 0
#             clade1 = clades[min_i]
#             clade2 = clades[min_j]
#             clade1.branch_length = dm[min_i, min_j] / 2.0
#             clade2.branch_length = dm[min_i, min_j] - clade1.branch_length
#             inner_clade = BaseTree.Clade(None, "Inner")
#             inner_clade.clades.append(clade1)
#             inner_clade.clades.append(clade2)
#             clades[0] = inner_clade
#             root = clades[0]
# 
#             return BaseTree.Tree(root, rooted=False)
#         while len(dm) > 2:
#             print("Joining pop#", len(dm))
#             # calculate nodeDist
#             for i in range(0, len(dm)):
#                 node_dist[i] = 0
#                 for j in range(0, len(dm)):
#                     node_dist[i] += dm[i, j]
#                 node_dist[i] = node_dist[i] / (len(dm) - 2)
# 
#             # find minimum distance pair
#             min_dist = dm[1, 0] - node_dist[1] - node_dist[0]
#             min_i = 0
#             min_j = 1
#             for i in range(1, len(dm)):
#                 for j in range(0, i):
#                     temp = dm[i, j] - node_dist[i] - node_dist[j]
#                     if min_dist > temp:
#                         min_dist = temp
#                         min_i = i
#                         min_j = j
#             # create clade
#             clade1 = clades[min_i]
#             clade2 = clades[min_j]
#             inner_count += 1
#             inner_clade = BaseTree.Clade(None, "Inner" + str(inner_count))
#             inner_clade.clades.append(clade1)
#             inner_clade.clades.append(clade2)
#             # assign branch length
#             clade1.branch_length = (
#                 dm[min_i, min_j] + node_dist[min_i] - node_dist[min_j]
#             ) / 2.0
#             clade2.branch_length = dm[min_i, min_j] - clade1.branch_length
# 
#             # update node list
#             clades[min_j] = inner_clade
#             del clades[min_i]
# 
#             # rebuild distance matrix,
#             # set the distances of new node at the index of min_j
#             for k in range(0, len(dm)):
#                 if k != min_i and k != min_j:
#                     dm[min_j, k] = (
#                         dm[min_i, k] + dm[min_j, k] - dm[min_i, min_j]
#                     ) / 2.0
# 
#             dm.names[min_j] = "Inner" + str(inner_count)
#             del dm[min_i]
#             # dm.matrix = withoutIndices(dm.matrix, (min_i))
# 
#         # set the last clade as one of the child of the inner_clade
#         root = None
#         if clades[0] == inner_clade:
#             clades[0].branch_length = 0
#             clades[1].branch_length = dm[1, 0]
#             clades[0].clades.append(clades[1])
#             root = clades[0]
#         else:
#             clades[0].branch_length = dm[1, 0]
#             clades[1].branch_length = 0
#             clades[1].clades.append(clades[0])
#             root = clades[1]
# 
#         return BaseTree.Tree(root, rooted=False)
# =============================================================================
    
    def alleleFrequency(self):
        """Print an allele frequency table."""
        for pop in self.excelData:
            for locus in pop.loci:
                for allele, freq in pop.frequency(locus).items():
                    print(pop.name, locus, allele, round(freq, ndigits=4))

    def executeCommand(self):
        """Excute the workflow, display on command console."""
        try:
            # Load an excel file contain VNTR data
            self.loadExcelData()
            query = input('\nIs the displayed information correct? [y/n] ')
            if query.lower() != 'y':
                raise CancelException(
                    "VNTR information considered incorrect by user.")
            else:
                # Build a phylogenetic tree
                self.buildTree(formula='Cavalli')
                print('\nNeighbor-Joining tree constructed.')
                # Save the tree in specified file
                # self.saveTreeFile() # TODO
                print(self.tree)
                input('\nPress Enter to exit.')
        except CancelException as e:
            print(f'\n***{e.message}***')
            input('Press Enter to exit.')
    
    def executeGUI(self):
        """Excute the workflow, display on GUI."""
        pass

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
        self.loci = {}  # dict {locus:{allele:count,},}

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



class Tree():
    parent = None
    
    def __init__(self, left = None, right = None, dleft = 1, dright = 1):
        self.left  = left
        left.parent = self
        self.right = right
        left.dparent = dleft
        right.parent = self
        self.dleft  = dleft
        right.dparent = dright
        self.dright = dright
        

    
    def __repr__(self):
        head = lambda arr: arr[0]
        tail = lambda arr: arr[1:]
        s = "┬" + ("─" * int(self.dleft))
        if self.left:
            left = repr(self.left).splitlines()
            s += "" + head(left) + "\n"
            for line in tail(left):
                s += "│"  + (" " * int(self.dleft))
                s += line + "\n"
        if self.right:
            right = repr(self.right).splitlines()
            s += "└" + ("─" * int(self.dright))
            s += head(right) + "\n"
            for line in tail(right):
                s += " "  + (" " * int(self.dright))
                s += line + "\n"
        return s
    
    def toClades(self):
        
        left = None
        right = None
        if self.left:
            left = self.left.toClades()
        if self.right:
            right = self.right.toClades()
        try:
            return BaseTree.Clade(self.dparent, None, [left, right])
        except AttributeError:
            return BaseTree.Clade(1, None, [left, right])
        root = None
        return BaseTree.Tree(root, rooted=False)
    
    def __eq__(self, other):
        if other.__class__.__name__ != "Tree":
            return False
        if self.left == other.left:
            return  self.left   == other.left   and \
                    self.right  == other.right  and \
                    self.dleft  == other.dleft  and \
                    self.dright == other.dright
        else:
            return  self.left   == other.right  and \
                    self.right  == other.left   and \
                    self.dleft  == other.dright and \
                    self.dright == other.dleft
    def contains(self, other):
        if id(self) == id(other):
            return True
        else:
            if other.parent == None:
                return False
            else:
                return self.contains(other.parent)
    def distanceTo(self, other):
        if id(self) == id(other):
            return 0
        elif self.contains(other):
            return other.dparent + self.distanceTo(other.parent)
        elif other.contains(self):
            return self.dparent + other.distanceTo(self.parent)
        elif self.parent and other.parent:
            return self.dparent + other.dparent + self.parent.distanceTo(other.parent)
        else:
            raise LookupError("Nodes are not in the same tree!")

class Leaf(Tree):
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return self.name
    def __eq__(self, other):
        return self.name == other.name
    def toClades(self):
        return BaseTree.Clade(self.dparent, self.name)

if __name__ == '__main__':
    
    njtree = NJTreeConstructor()
    njtree.executeCommand()
    