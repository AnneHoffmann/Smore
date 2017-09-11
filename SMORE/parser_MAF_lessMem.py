#this program takes multiple sequence alignments in MAF format and parses them.
#It then makes BED files for each species that contain the allignments for that species.
import subprocess
import datetime #module that can give the current day month and year
from os import listdir
from os.path import isfile, isdir, join
from bisect import bisect_right #binary search
from dataStructures_lessMem import Block, Chromosome
from parser_infernal import parseInfernal

'''
Order of the data structures:

   Chromosomes -> Blocks

   Each chromosome has a list of the multiple sequence allignments that exist on that chromosome

Definitions:
   Reference species: The species that the multiple sequence allignment was done against
   Reference sequence: The block objects that are from the reference species, contain a list of allignment sequences
   Allignment sequence: The block objects that are from all the other non-reference species. Dont know what their reference seq is

Procedure:
   Read the MAF file line by line. Write each allignnet to a temp file corresponding to it's species.
   Once all MAF files have been read, sort the temp files so the allignments are grouped by chromosome. Read temp files one
   chromosome at a time. As each allignment is read it is sorted via binary search into it's correct place on the chromosome
   object's list. If the chromosome is from the reference species we check if the allignments overlap. Once all of the sequences have been added to the chromosome's 
   list we search the reference species for overlapping reference sequences and mark all the sequences alligned to them
   as being overlapping.
 
   We only look for overlapping blocks in the reference species.
   Once that is done we remove all overlapping reference blocks / sequeces  and the sequences from non-reference
   species that are alligned with them.
'''
def catchExceptions(fileName):#, outputDir):
    '''
    checks if file is in correct format( '.maf', '.maf.gz', '.maf.Z', '.maf.bz2')
    checks if file is zipped
    checks if file can be opened
    '''

    #proper extension
    if fileName.endswith(('.maf.gz','.maf.Z','.maf.bz2','.maf')):
        pass
    else:
        raise Exception("{} not of maf format. Use a .maf file or gziped .maf file('.maf.gz','.maf.Z','.maf.bz2')".format(fileName))

    #can file be opened
    try:
        f = open(fileName, 'r')
    except IOError:
        raise Exception("cannot open {}".format(fileName))
    else:
        f.close()
        return


def maf2TempWrapper(mafDir, outputDir, repoDir, listOfSpecies):
    '''
    outputDir = ../test/Insects/Output/
    Assumptions:
        The first species in the list is the reference species

    open a file for each species given. give those files to the maf to temp parser
    close all the files
    '''

    #print('starting temp')
    if 'temp' in listdir(outputDir):
        subprocess.call("rm -r "+outputDir+'temp', shell=True)
    subprocess.call("mkdir "+outputDir+'temp', shell = True)
    tempOutputDir = outputDir + 'temp/'

    #Make temp files
    speciesFiles = list()
    for species in listOfSpecies:
        tempBedFile = open(tempOutputDir+species+'_temp.bed', 'w')
        tempBedFile.close()
        speciesFiles.append(tempBedFile.name)

    #NEW PIPING OPERATION

    blockNum = 0
    iteration = 0

    mafFileNames = [join(mafDir, f) for f in listdir(mafDir) if isfile(join(mafDir,f))]

    #cwd = os.getcwd()
    #print("cwd: {}".format(cwd))
    #if cwd.endswith(
    for name in mafFileNames:
        catchExceptions(name)
        
        #use Popen to catch stdOutput of reading file, stdOutput being the current block number
        proc = subprocess.Popen("zcat "+name+" | python3 "+repoDir+"maf2TempBed.py "+str(blockNum)+" "+tempOutputDir, shell=True, universal_newlines=True, stdout=subprocess.PIPE)

        #print('read: {}'.format(name))
        stdOutput, errors = proc.communicate()
        #print('output: {}'.format(stdOutput))
        try:
            oldblockNum = blockNum
            blockNum = int(stdOutput)
        except ValueError:
            blockNum = oldblockNum
    #End of piping operation
    
    for _file in speciesFiles:
        #sort -k1,1 -k2,2n outputDir/dm6_temp.bed > outputDir/dm6_temp_sorted.bed
        subprocess.call("sort -k1,1 "+_file+" > "+_file.replace('.bed','_sorted.bed', 1), shell=True)
        #print('sorted: {}'.format(_file))
        
    #returnList is a list of file paths.
    #Initialized with the reference species as the first element
    #however this is done with literals so be careful if changing the file name schema
    sortedTempList = [join(tempOutputDir, listOfSpecies[0]+'_temp_sorted.bed')]
    for _file in listdir(tempOutputDir):
        fileName = join(tempOutputDir, _file)
        if fileName not in sortedTempList:
            if fileName.endswith("_sorted.bed"):
                sortedTempList.append(fileName)
            else:
                subprocess.call("rm "+fileName, shell=True)

    #print('ending temp')
    return sortedTempList

def parseTempWrapper(tempFiles, listOfListOfGenes, outputDir, threshold, infernalVersion):
    '''
    Assumptions:
        The first file in tempFiles must be the file for the reference Species. As a result this means
        the first species in the listOfSpecies given to the maf parser must be the reference species.

        The tempFiles must be sorted by chromosome
    '''

    #print('start parsing temp')
    bedFileDir = outputDir+'bed/'
    if 'bed' in listdir(outputDir):
        subprocess.call("rm -r "+bedFileDir, shell=True)
    subprocess.call("mkdir "+bedFileDir, shell=True)

    geneFileDir = outputDir+'genes/'
    if 'genes' in listdir(outputDir):
        subprocess.call("rm -r "+geneFileDir, shell=True)
    subprocess.call("mkdir "+geneFileDir, shell=True)        
    
    #print("parsing temp files...")    
    overlappingBlockNums = set()

    now = datetime.datetime.now()
    
    for tempBlockFileName in tempFiles:

        tempBlockFile = open(tempBlockFileName, 'r')
        
        species = tempBlockFile.name.split('/')[-1].split('_')[0]
#        print("species: {}".format(species))
#        print("len geneList: {}".format(len(listOfListOfGenes)))

        finalBlockFileName = bedFileDir+species+'.bed'
        finalBlockFile = open(finalBlockFileName, 'w')
        start = finalBlockFile.write("# Created on "+str(now.day)+"/"+str(now.month)+"/"+str(now.year)+" (day/month/year)\n"+\
                                     "\n")
        finalBlockFile.seek(start)
        geneFile = open(geneFileDir+species+'.bed', 'w')
        start = geneFile.write("# Created on "+str(now.day)+"/"+str(now.month)+"/"+str(now.year)+" (day/month/year)\n"+\
                               "# "+infernalVersion+"\n"+\
                               "\n")
        geneFile.seek(start)

        listOfGenes = None
        for _list in listOfListOfGenes:  
            if _list[0].species == species:
                listOfGenes = _list

        overlappingBlockNums = parseTemp(tempBlockFile, finalBlockFile, geneFile, listOfGenes, overlappingBlockNums, threshold)
        #print('done')

        tempBlockFile.close()
        finalBlockFile.close()
        geneFile.close()

        subprocess.call('gzip '+finalBlockFileName, shell=True)
        subprocess.call('gzip '+tempBlockFileName, shell=True)
        #subprocess.call('rm '+tempBlockFileName, shell=True)
    '''
    if len(listOfGenes) > 0:
        listOfSC = []
        for i in range(len(listOfGenes)):
            if (listOfGenes[i].species, listOfGenes[i].chromosome) not in listOfSC:
                listOfSC.append((listOfGenes[i].species, listOfGenes[i].chromosome))
        print('unsorted genes')
        for tup in listOfSC:
            if tup[0] == 'droMoj3':
                print(tup)

    else:
        print('no unsorted genes')
    '''

def parseTemp(tempFile, finalFile, geneFile, listOfGenes, overlapSet, threshold):
    '''
    Read a single chromosome, check all the overlaps and gene collisions, write the chromosome to
    the finalFile (output file). Read the next chromosome on the temp file, continue as such till all
    of the chromosomes are read. If the 
    '''
    #print("len listOfGenes begin parseTemp: {}".format(len(listOfGenes)))
    #print("threshold parsetemp: {}".format(threshold))
    chromo = None
    #print("parsing temp file...")
    for line in tempFile:
        lineCont = line.split()
        species_blockNum = lineCont[1].split('_')
#        print("species blockNum read Chromo: {}".format(species_blockNum))
        
        if chromo == None:
            chromo = Chromosome(lineCont[0], species_blockNum[0])
            if lineCont[5] == 'True':
                chromo.isReference = True

 #       print("chromo.name: {}".format(chromo.name))
 #       print("lineCont[0]: {}".format(lineCont[0]))
        #if we are still on the same chromosome
        if lineCont[0] == chromo.name:
        #                        start pos       length           strand       blockNum                  score
            chromo.add(Block(int(lineCont[2]),int(lineCont[3]),lineCont[4],int(species_blockNum[1]), int(lineCont[6])))
            if threshold > 0:
                chromo.listOfScores.append(int(lineCont[6]))

        #if we reached the end of the chromosome.
        #check geneoverlaps, write genes that are from the last chromosome.
        #remove overlaps from the last chromosome, write it to the bed file.
        #overwrite 'chromo' to the new chromosome object.
        else:
#            print("len chromoListOfScores bef readChromo: {}".format(len(chromo.listOfScores)))
#            print("len listOfGenes bef readChromo: {}".format(len(listOfGenes)))
            overlapSet, listOfGenes = readChromo(finalFile, geneFile, listOfGenes, overlapSet, threshold, chromo)
            #overwrite chromo to new chromosome
            chromo = Chromosome(lineCont[0], species_blockNum[0])

    #read the last chomo in the file
    overlapSet, listOfGenes = readChromo(finalFile, geneFile, listOfGenes, overlapSet, threshold, chromo)

    listOfGenes = listOfGenes + listOfGenes
    #print("len geneList parseTemp: {}".format(len(listOfGenes4)))
    
    if len(listOfGenes) > 0:
#        print('writing orphin')
#        print("species: {}".format(chromo.species))
#        print('chromos:  {}'.format([gene.chromosome for gene in listOfGenes4]))
        for gene in listOfGenes:
            fivePrime, threePrime = None, None
            if not gene.ownGene:
                geneFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
                               gene.chromosome, chromo.species+"_"+str(gene.blockNum), gene.s,\
                               gene.getEndPos(), gene.strand, fivePrime, threePrime, gene.structure, gene.sequence, gene.score))
            else:
                geneFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
                               gene.chromosome, chromo.species+"_"+str(gene.blockNum), gene.s,\
                               gene.getEndPos(), gene.strand, fivePrime, threePrime, gene.structure, gene.sequence, \
                               gene._type, gene.pseudo, gene.comment))

    return overlapSet

def readChromo(finalFile, geneFile, listOfGenes, overlapSet, threshold, chromo):

#    print("len geneList readChromo: {}".format(len(listOfGenes)))
#    print("len chromo.listOfScores readchromo: {}".format(len(chromo.listOfScores)))
    
    if len(chromo.listOfScores) > 0:
        chromo.listOfScores.sort()
        thresholdIndex = int(len(chromo.listOfScores)/100*threshold)
        if thresholdIndex < 0:
            thresholdIndex = 0
        elif thresholdIndex >= len(chromo.listOfScores) and len(chromo.listOfScores) > 0:
            thresholdIndex = len(chromo.listOfScores) -1
            thresholdScore = chromo.listOfScores[thresholdIndex]
    else:
        thresholdScore = float('-inf')#if threshold is 0 set thresholdScore to - infinity
                
            #We check the entire gene list for every chromosome for every species but
            #we make a new list of the unwritten genes and use that for the next iteration
            #so the list always gets shorter with each iteration. And we can see if any genes are not sorted
            
    remainingGenes = []
            
    for i in range(len(listOfGenes)):
        gene = listOfGenes[i]
        if gene.chromosome == chromo.name:
            chromo.checkGene(gene)
            fivePrime, threePrime = chromo.getAdjBlock(gene)
            if not gene.ownGene:
                geneFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
                               chromo.name, chromo.species+"_"+str(gene.blockNum), gene.s,\
                               gene.getEndPos(), gene.strand, fivePrime, threePrime, gene.structure, gene.sequence, gene.score))
            else:
                geneFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
                               gene.chromosome, chromo.species+"_"+str(gene.blockNum), gene.s,\
                               gene.getEndPos(), gene.strand, fivePrime, threePrime, gene.structure, gene.sequence, \
                               gene._type, gene.pseudo, gene.comment))
                            
        else:
            remainingGenes.append(gene)
            
    listOfGenes = remainingGenes
    remainingGenes = None #free up memory by not having two copies of the same list

            #replaced deleting with creating new list. Remove had O(n) cost + n operations
            #making new list means appending O(1) + n operations
    if chromo.isReference:
        writeList = []
        for block in chromo.listOfMultiZ:
            if block.Overlap:
               overlapSet.add(block.blockNum)
            elif block.score < thresholdScore:
               continue
            else:
               #only add blocks that are above threshold and dont overlap
               writeList.append(block)
    else:
                #only blocks above threshold, have valid reference block and dont overlap gene are added
        writeList = [block for block in chromo.listOfMultiZ if not block.Overlap and block.blockNum not in overlapSet and block.score > thresholdScore]

    #free up memory in chromo
    chromo.listOfMultiZ = None

    #writing blocks
    for i in range(len(writeList)):
                
        block = writeList[i]
        if i-1 < 0:
            fivePrime = None
        else:
            fivePrime = writeList[i-1].blockNum
        if i+1 >= len(writeList):
            threePrime = None
        else:
            threePrime = writeList[i+1].blockNum
        finalFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
                        chromo.name, chromo.species+"_"+str(block.blockNum), block.s,\
                        block.getEndPos(), block.strand, fivePrime, threePrime))

    return overlapSet, listOfGenes

def maf2bed(mafDir, outputDir, repoDir, geneObjs, listOfSpecies, threshold, infernalVersion):
    '''
    Assumptions:
        The first string in listOfSpecies must correspond to the reference species.
    '''
    tempFiles = maf2TempWrapper(mafDir, outputDir, repoDir, listOfSpecies)

    parseTempWrapper(tempFiles, geneObjs, outputDir, threshold, infernalVersion)

'''    
#if this python file is called from the command line
if __name__ == "__main__":

    
    listOfSpecies = ['dm6', 'anoGam1', 'apiMel4', 'droAlb1', 'droAna3', 'droBia2', 'droBip2', 'droEle2', 'droEre2', 'droEug2',\
                     'droFic2', 'droGri2', 'droKik2', 'droMir2', 'droMoj3', 'droPer1', 'droPse3', 'droRho2', 'droSec1', 'droSim1',\
                     'droSuz1', 'droTak2', 'droVir3', 'droWil2', 'droYak3', 'musDom2', 'triCas2'] #------------------------Testing -----------------------------

    maf2TempWrapper("/scr/rum/cameron/maf","/homes/biertruck/cameron/Desktop/Project_May_June_2016/test/Insects/Output8/",listOfSpecies)
'''
