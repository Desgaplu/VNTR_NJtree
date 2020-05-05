# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 22:39:15 2020

@author: pldesgagne

Description:
    Blah # TODO
"""


class NJTreeConstructor():
    """
    # TODO
    """
    def __init__(self):

        self.excelData = None
        self.distanceMatrix = None
        self.tree = None
        
    
    def loadExcelData(self, source:str):
        """
        Load the data from excel sheet into a list of Pop()
        which contains name, locus, alleles for each samples.
        Saves results into class variable self.excelData
        
        Input:
            source: file name and address for the file containing the samples
                name and their allele value for each locus.
                format in each colomn should be:
                    names locus values # headers
                    name locus1 value1
                    name locus2 value2
                    
        Returns
            None.
        """
        # try reading file if it exist
        # load data into list of pop then into self.excelData
        # scan the data
        
        pass
    
    def __scanExcelData(self, data):
        """
        Parse the excel data and show different information which may help the
        user find potential error in the provided excel files such as typos in
        names or wrong sample number.
        
        Input:
            data: the data to be scanned in format [Pop()]
            
        Returns
            None.
        """
        # find inconsistencies and print a warning (compare loci names)
        # print a data feed back (number of samples, loci etc):
        # number of samples loaded.
        # list of loci names and their relative count
        # max number of  of a single loci
            # locus names are lowered to prevent cap typos
        # pop with less than 4 loci (potential typo in pop name)
        pass # TODO
    
    def saveTreeFile(self, destination:str):
        """
        Save the tree data into self.destination file .___ # TODO
        
        Input:
            destination: file name and address for the tree file created by
                this script.
                
        Returns
            None.
        """
        pass # TODO
    
    def buildTree(self, data=None):
        """
        Takes the VNTR data and return a Neighbor-Joining tree.
        Uses the Cavalli-Sforza chord distance formula to create a distance
        matrix then the Biopython Phylo package to build the neighbor tree
        
        Input:
            data (optionnal): the VNTR data to be analysed. If none is passed,
                the self.excelData from the instance is used.
        
        Returns
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
        
        Input:
            data: the VNTR data loaded by the loadExcelData method.
            
        Return:
            DistanceMatrix from the Biopython Phylo TreeConstruction module
        """
        pass # TODO
    
    def __neighbor(self, matrix):
        """
        Apply the neighbor joining algo on a distance matrix and returns 
        a tree

        Input:
            matrix : DistanceMatrix from the Biopython Phylo TreeConstruction 
                module

        Returns
            Tree from the Biopython Phylo BaseTree module
        """
        pass # TODO
    
    
class Pop():
    """
    Representation of a single sample.
    Attributes:
            name: Name of sample in string format
            loci: a dict of each locus name containing alleles values.
                {loci_name:{allele_value:count,},}
    """
    def __init__(self, name):
        """
        Input:
            name: Name of sample in string format.
        """
        self.name = name
        self.loci = {} # dict {name:{allele:count,},}
        
    def addLocus(self, locus:str, allele:int):
        """
        Adds the locus and the allele value in the loci dict.
        If the allele is already present for that loci, +1 to count.
        """
        if locus not in self.loci:
            self.loci[locus] = {allele:1}
        elif allele in self.loci[locus]:
            self.loci[locus][allele] += 1
        else:
            self.loci[locus][allele] = 1
    
    def frequency(self, locus):
        """
        Return a dict {allele:frequency} of the frequency of each allele
        for the specified locus.
        
        Ex.: for 3 alleles {650:2, 720:1} => {650:0.6666, 720:0.3333} 
        (Divides the count of each allele by the total number of allele for
        the specified loci.)
        """
        nb_allele = 0
        for allele in self.loci[locus]:
            nb_allele += self.loci[locus][allele]
        return {i:(self.loci[locus][i]/nb_allele) for i in self.loci[locus]}

def testUnit():
    sample = Pop('R1105')
    sample.addLocus('Tr1', 650)
    sample.addLocus('Tr1', 650)
    sample.addLocus('Tr1', 720)
    sample.addLocus('Tr2', 150)
    # assert test_bool, "Test failed comment"
    # print("All test passed")

# =============================================================================
# if __name__ == '__main__':
#     
#     njtree = NJTreeConstructor()
#     source = input("Enter the excel file name: ")
#     njtree.loadExcelData(source)
#     query = input('Is the displayed information correct? [y/n] ')
#     if query.lower() != 'y':
#         print('\nCancelling...')
#         input('Press Enter to exit.')
#     else:
#         njtree.buildTree()
#         print('Neighbor-Joining tree constructed.')
#         destination = input('Enter a name for the tree file: ')
#         njtree.saveTreeFile(destination)
#         print('Tree saved!')
#         input('Press Enter to exit.')
# =============================================================================
