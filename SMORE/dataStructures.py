
from bisect import bisect_right #binary search returns index where  (O(n) because python insertion is O(n))


class Sequence:
    '''
    Species:    A string containing the name of the species
    Chromosome: A string containing the name of the chromosome
    StartPos:   An int of the nucleotide position in the chromosome sequence
    length:     An int of the length of the sequence
    strand:     A + or - char depending if the alignment is on the positive or negative strand.
                  If the strand is + add the len to start to get the end pos. If the strand is -
                  subtract the len from the start to get the end.
    blockNum:   The alignment block number across all the multiple sequence alignments
    Overlap:    A bool, True if sequence overlaps with another sequence, False if not
    geneOverlap:A bool, True if sequence overlaps with a gene, False if not
    fivePrime:  The block number of the upstream block
    threePrime: The block number of the downstream block
    refSeq:     A bool, True if this sequence is the reference sequence, False if just an alligned sequence.
    isGene:     A bool, True if this sequence is a gene, False if its an alignment
    '''

    def __init__(self, species, chromosome, startPos, length, strand, refSeq, isGene, blockNum = None):
        self.species = species
        self.chromosome = chromosome
        self.s = startPos
        self.length = length
        self.strand = strand
        self.blockNum = blockNum
        self.Overlap = False
        self.geneOverlap = False
        self.fivePrime = None
        self.threePrime = None
        self.refSeq = refSeq
        if self.refSeq:
            self.listOfAllignedSeqs = list()
        self.isGene = isGene
        
        #self.sequence = None #maybe for later when we check uniquness

    def __lt__(self, multiZ):
        '''
        called when python evaluates Sequence < Sequence

        allows the blocks to be arranged in order of their position 
        relitive to the 5' end of the + strand.
        '''
        if self.strand == '+':
            comparNum1 = self.s
        else:
            comparNum1 = self.getEndPos()
        if multiZ.strand == '+':
            comparNum2 = multiZ.s
        else:
            comparNum2 = multiZ.getEndPos()

        if comparNum1 < comparNum2:
            return True
        return False

    def __gt__(self, multiZ):
        '''
        called when python evaluates Sequence > Sequence
        '''
        if self.strand == '+':
            comparNum1 = self.s
        else:
            comparNum1 = self.getEndPos()
        if multiZ.strand == '+':
            comparNum2 = multiZ.s
        else:
            comparNum2 = multiZ.getEndPos()

        if comparNum1 > comparNum2:
            return True
        return False

            
    def getEndPos(self):
        '''
        returns the end position of the allignment

        if we are on the + strand we read left to right so add. If we are on the - strand
        we read right to left so we subtract.
        '''
        if self.strand == "+":
            return self.s + self.length
        return self.s - self.length

        
class Allignments:
    '''
    container to hold all of the species objects in an allignment
    '''
    def __init__(self):
        self.listOfSpecies = list()

    def add(self, multiZ):
        for species in self.listOfSpecies:
            if multiZ.species == species.name:
                species.add(multiZ)
                return
        else:
            new_species = Species(multiZ.species)
            new_species.add(multiZ)
            self.listOfSpecies.append(new_species)

            
    def checkGene(self, gene):
        for species in self.listOfSpecies:
            if gene.species == species.name:
                species.checkGene(gene)
                return
        else:
            raise Exception("gene does not have properly formatted species name \n"+
                            "name given was {}. Acceptable names are {}".format(gene.species, *[species.name for species in self.listOfSpecies]))

    def getAdjBlocks(self, gene):
        for species in self.listOfSpecies:
            if gene.species == species.name:
                #print("species gene is in: {}".format(species.name))
                nums = species.getAdjBlocks(gene)
                #print("species nums: {}".format(nums))
                return nums
        
        else:
            raise Exception("gene does not have properly formatted species name \n"+
                            "name given was {}. Acceptable names are {}".format(gene.species, *[species.name for species in self.listOfSpecies]))

class Species:
    '''
    named container that holds the chromosome objects of a given species
    '''

    def __init__(self, name):
        self.name = name
        self.listOfChromosomes = list()

    def add(self, multiZ):
        for chromosome in self.listOfChromosomes:
            if multiZ.chromosome == chromosome.name:
                chromosome.add(multiZ)
                return
                
        else:
            new_chromo = Chromosome(multiZ.chromosome)
            new_chromo.add(multiZ)
            self.listOfChromosomes.append(new_chromo)

    def checkGene(self, gene):
        for chromosome in self.listOfChromosomes:
            if gene.chromosome == chromosome.name:
                chromosome.checkGene(gene)
                return
        else:
            pass #---------------------------------------------------FOR TESTING PURPOSES ONLY ---------------------------------------------------------
            #raise Exception("gene does not have properly formatted chromosome name \n"+
             #               "name given was {}. Acceptable names are {}".format(gene.chromosome, *[chrom.name for chrom in self.listOfChromosomes]))

    def getAdjBlocks(self, gene):
        for chromosome in self.listOfChromosomes:
            if gene.chromosome == chromosome.name:
                #print("chromosome gene is in: {}".format(chromosome.name))
                nums = chromosome.getAdjBlock(gene)
                #print("chromosome level: {}".format(nums))
                return nums
        
        else:
            #print("gene not in any chromosome we know of")
            return ("None", "None")
            #pass ---------------------------------------------------FOR TESTING PURPOSES ONLY ---------------------------------------------------------
            #raise Exception("gene does not have properly formatted chromosome name \n"+
             #               "name given was {}. Acceptable names are {}".format(gene.chromosome, *[chrom.name for chrom in self.listOfChromosomes]))

            
class Chromosome:

    def __init__(self, name):
        self.name = name
        self.listOfMultiZ = list()

    def add(self, multiZ):
        '''
        O(n) running time
        '''
        i = bisect_right(self.listOfMultiZ, multiZ)
        self.listOfMultiZ.insert(i, multiZ)
        if multiZ.refSeq == True:
            self.checkOverlap(i, multiZ)

    def checkOverlap(self, i, multiZ):

        self.leftCheck(i, multiZ)
        self.rightCheck(i, multiZ)
        return

    def checkGene(self,gene):
        '''
        almost the same as add. Does not insert gene into the list of Multiz
        '''
        
        if gene.isGene == False:
            raise Exception("Sequence given was not a gene.")
        i = bisect_right(self.listOfMultiZ, gene)
        self.checkOverlapGene(i, gene)

    def checkOverlapGene(self, i, gene):
        
        self.leftCheck(i, gene)
        self.rightCheck(i-1, gene)#since we arent inserting the gene we need to check the block
                                  #that the gene would displace
        return

    def getAdjBlock(self, gene):
        if gene.isGene == False:
            raise Exception("Sequence given was not a gene.")
        i = bisect_right(self.listOfMultiZ, gene)
        #print("index in chromosome list: {}".format(i))
        if i-1 >= 0:
            #print("middle of chromosome")
            #print("left object: {}".format(self.listOfMultiZ[i-1]))
            #print("obj blockNum: {}".format(self.listOfMultiZ[i-1].blockNum))
            #print("obj fivePrime: {}".format(self.listOfMultiZ[i-1].fivePrime))
            adjBlockNums = (self.listOfMultiZ[i-1].blockNum, self.listOfMultiZ[i].blockNum)
            #print("return Value: {}".format(adjBlockNums))
        else:
            #print("end of chromosome")
            adjBlockNums = (None, self.listOfMultiZ[i].blockNum)
        return adjBlockNums
    
#Functions for checking the neighbours
        
            
    def leftCheck(self, i, multiZ):
        #is the block to the immediate left overlapping
        #the block can have the first basepair in common with the immediate left's last

        if i -1 >= 0:#return immediately if block is on the left end of the chromosome
            left = self.listOfMultiZ[i-1]

            #conditions for overlapping
            if multiZ.strand == '+':
                    
                if left.strand == '+' and left.getEndPos() > multiZ.s:
                    parameterSet(left, multiZ)
                elif left.strand == '-' and left.s > multiZ.s:
                    parameterSet(left, multiZ)
                
            else:#multiZ is a - strand
                if left.strand == '+' and left.getEndPos() > multiZ.getEndPos():
                    parameterSet(left, multiZ)
                elif left.strand == '-' and left.s > multiZ.getEndPos():
                    parameterSet(left, multiZ)

        return

    def rightCheck(self, i, multiZ):
        #iterativly search of all the blocks to the right.
        #Searches until a non-over lapping block is found.
        if i +1 < len(self.listOfMultiZ):#return immediately if block is on the right end of the chromsome
            k = i
        
            while True:
                k+=1
                right = self.listOfMultiZ[k]

                #conditions for overlapping
                if multiZ.strand == '+':
                    if right.strand == '+' and right.s < multiZ.getEndPos():
                        parameterSet(right, multiZ)
                        continue#return to top of loop
                    elif right.strand == '-' and right.getEndPos() < multiZ.getEndPos():
                        parameterSet(right, multiZ)
                        continue
                
                else:#multiZ is a - strand
                    if right.strand == '+' and right.s < multiZ.s:
                        parameterSet(right, multiZ)
                        continue
                    elif right.strand == '-' and right.getEndPos() < multiZ.s:
                        parameterSet(right, multiZ)
                        continue        
                return
        return

def parameterSet(overlapObj, multiZ):
    if multiZ.isGene == True:
        overlapObj.geneOverlap = True
#       multiZ.geneOverlap = True
    else:
        overlapObj.Overlap = True
        multiZ.Oberlap = True

