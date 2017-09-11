
from bisect import bisect_right #binary search returns index where  (O(n) because python insertion is O(n))



class Block:
    '''
    StartPos:   An int of the nucleotide position in the chromosome sequence
    length:     An int of the length of the sequence
    strand:     A + or - char depending if the alignment is on the positive or negative strand.
                  If the strand is + add the len to start to get the end pos. If the strand is -
                  subtract the len from the start to get the end.
    blockNum:   The alignment block number across all the multiple sequence alignments
    Overlap:    A bool, True if block overlaps with another block. False if not
    		Only reference sequence blocks are checked.
    score:      An int read from the maf files. Used when the -q flag is given for determining
		Which blocks should be kept and which should be thrown away.
    '''

    def __init__(self, startPos, length, strand, blockNum, score):

        self.s = startPos
        self.length = length
        self.strand = strand
        self.blockNum = blockNum
        self.score = score
        self.Overlap = False

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

class Gene (Block):
    '''
    Species:      A string containing the name of the species
    Chromosome:   A string containing the name of the chromosome
    StartPos:     An int of the nucleotide position in the chromosome sequence
    length:       An int of the length of the sequence
    strand:       A + or - char depending if the alignment is on the positive or negative strand.
                    If the strand is + add the len to start to get the end pos. If the strand is -
                    subtract the len from the start to get the end.
    geneNum:      The unique int asociated with each gene
    structure:    A string of the ascii representation of the secondary structure of the genes found by cmsearch
    sequence:     A string of the DNA sequence of the gene in single letter bases (ex: TGCTA)

    if searching for genes:
        Score:    An int of the Infernal score for the element. Used for determingin pseudogenes later on. 
    if own genes:
        _type:    A string of the type of gene given (ex: 'tRNA' or 'heme_isoform').
        pseudo:   A bool if the gene given is a pseudogene or not.
        comment:  A string a comment the user wants to leave about a particular element.
    '''

    def __init__(self, species, chromosome, startPos, length, strand, geneNum, structure, sequence, score=None, ownGene=False, _type=None, pseudoGene=None, comment=None):
        self.species = species
        self.chromosome = chromosome
        self.s = startPos
        self.length = length
        self.strand = strand
        self.blockNum = geneNum
        self.sequence = sequence
        self.structure = structure
        self.ownGene = ownGene
        if not ownGene:
            self.score = score
        if ownGene:
            self._type = _type
            self.pseudo = pseudoGene
            self.comment = comment
        
    def __lt__(self, multiZ):
        '''
        called when python evaluates Sequence < Sequence

        allows the genes to be arranged in order of their position
        relitive to the 5' end of the + strand.
        '''

        return Block.__lt__(self, multiZ)


    def __gt__(self, multiZ):
        '''
        called when python evaluates Sequence > Sequence
        '''
        return Block.__gt__(self, mlultiZ)

    def getEndPos(self):
        '''
        returns the end position of the allignment
        '''
        return Block.getEndPos(self)

class Chromosome:
    '''
    name:         A string of the chromosome's name as it appears in the multiple sequence alignment
		    Ex: chr3
    species:      A string of the chromosome's species' name as it appears in the multiple sequence alingment
		    Ex: panTro
    isReference:  A bool. If this chromosome's species is the reference species True otherwise False
    listOfMultiZ: A list of all of the blocks in this chromosome in order of their positions. Where the
		  Beginning of the list is the block closest to the 5' end of the + strand
    listOfScores: A list of all of the block's scores from the multiple sequence alignment.
    '''

    def __init__(self, name, species):
        self.name = name
        self.species = species
        self.isReference = False
        self.listOfMultiZ = list()
        self.listOfScores = list()

    def add(self, multiZ):
        '''
        O(n) running time
        where n = number of elements in listOfMultiZ
        '''
        i = bisect_right(self.listOfMultiZ, multiZ)
        self.listOfMultiZ.insert(i, multiZ)
        if self.isReference == True:
            self.checkOverlap(i, multiZ)

    def checkOverlap(self, i, multiZ):

        self.leftCheck(i, multiZ)
        self.rightCheck(i, multiZ)
        return

    def checkGene(self,gene):
        '''
        almost the same as add. Does not insert gene into the list of Multiz
        '''
        
        if not isinstance(gene, Gene):
            raise Exception("Sequence given was not a gene.")
        i = bisect_right(self.listOfMultiZ, gene)
        self.checkOverlapGene(i, gene)

    def checkOverlapGene(self, i, gene):
        
        self.leftCheck(i, gene)
        self.rightCheck(i-1, gene)#since we arent inserting the gene we need to check the block
                                  #that the gene would displace, if it was being inserted
        return

    def getAdjBlock(self, gene):
        if not isinstance(gene, Gene):
            raise Exception("Sequence given was not a gene.")
        i = bisect_right(self.listOfMultiZ, gene)
        adjBlockNums = (None, None)
        if i-1 >= 0:
            adjBlockNums = (self.listOfMultiZ[i-1].blockNum, adjBlockNums[1])
        if i < len(self.listOfMultiZ):
            adjBlockNums = (adjBlockNums[0], self.listOfMultiZ[i].blockNum)
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
        #if i +1 < len(self.listOfMultiZ):#return immediately if block is on the right end of the chromsome
            #k = i
        k = i
        while True:
            if k +1 < len(self.listOfMultiZ):
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

def parameterSet(overlapObj, multiZ):
    overlapObj.Overlap = True
    if not isinstance(multiZ, Gene):
        multiZ.Overlap = True
