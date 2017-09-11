
#this program takes multiple sequence alignments in MAF format and parses them.
#It then makes BED files for each species that contain the allignments for that species.
#All blocks that are overlapping will not be included in the BED file.
import subprocess
from os import listdir
from os.path import isfile, join
from bisect import bisect_right #binary search
from dataStructures import Sequence, Allignments, Species, Chromosome
from parser_infernal import parseInfernal

'''
Order of the data structures:

   Allignments -> Species -> Chromosomes -> Sequences

   Allignments is a list of all the species
   Each Species has a list of chromosomes in it
   Each chromosome in a species has the sequences that exist on that chromosome

Definitions:
   Reference species: The species that the multiple sequence allignment was done against
   Reference sequence: The sequence objects that are from the reference species, contain a list of allignment sequences
   Allignment sequence: The sequence objects that are from all the other non-reference species. Dont know what their reference seq is

Procedure:
   Read the MAF file to aquire sequences. Add the sequences to the allignment. Adding the sequence sorts it into it's correct 
   species, chormosome and place on the chromosome. It also checks if the reference sequences overlap. Once all of the sequences
   have been added we search the reference species for overlapping reference sequences and mark all the sequences alligned to them
   as being overlapping. 
   We only look for overlapping blocks in the reference species.
   Once that is done we remove all overlapping reference blocks / sequeces  and the sequences from non-reference
   species that are alligned with them.
'''
def catchExceptions(fileName, outputDir):
    '''
    checks if file is in correct format( '.maf', '.maf.gz', '.maf.Z', '.maf.bz2')
    checks if file is zipped
    checks if file can be opened
    checks if file can be read

    unzips file if zipped
    '''
    
    name = fileName.split('/')[-1]#fileName without pathway ex: test1.maf.gz
    
    #proper extension
    if fileName.endswith('.maf'):
        pass

    elif fileName.endswith(('.maf.gz','.maf.Z','.maf.bz2')):

        unzippedDir = outputDir+'unzippedMAF'
        if 'unzippedMAF' in listdir(outputDir):
            subprocess.call("rm -r "+unzippedDir)
        subprocess.call('mkdir '+unzippedDir)
        
        #unzippedFileName = ('/'.join(fileName.split('/')[:-1]))+'/'+name.split('.')[0]+'.maf'
        unzippedFileName = unzippedDir +'/'+ name.split('.')[0]+'.maf'
        print("unzipping maf file...")
        subprocess.call('/usr/bin/zcat '+fileName+' > '+unzippedFileName, shell=True)
        print("done")
        fileName = unzippedFileName
    else:
        raise Exception("{} not of maf format. Use a .maf file or gziped .maf file('.maf.gz','.maf.Z','.maf.bz2')".format(fileName))

    #can unzipped file be opened + read
    try:
        f = open(fileName, 'r')
        f.readline()
    except IOError:
        raise Exception("cannot open {}".format(fileName))
#    except UnicodeDecodeError:
#        raise Exception("File not of maf format. Use an unziped .maf file")
    else:
        return f

def parseMAF(fileNames, listOfGenes, outputDir, allignments=Allignments()):
    '''
    for output to be consistent files in the directory need to be in order
    '''
    
    refSpecies = None
    print("parsing MAF files...")
    for fileName in fileNames:
        
        print("unzipping if needed...")
        f = catchExceptions(fileName, outputDir)
        print("done")
        print("adding blocks in maf files..."
        
        lastLineA = False
        currentRefSeq = None
        for line in f:
            if line[0] == 'a':
                lastLineA = True
                currentRefSeq = None
            if line[0] == 's':
                buf = line.split()
                ID = buf[1].split('.')
                    
                #if the previous line is 'a' the next line is the reference genome
                if lastLineA == True:
                    currentRefSeq = Sequence(ID[0],ID[1],int(buf[2]),int(buf[3]),buf[4],lastLineA, False)
                    allignments.add(currentRefSeq)
                    if refSpecies == None:
                        refSpecies = ID[0]
                else:
                    new_seq = Sequence(ID[0],ID[1],int(buf[2]),int(buf[3]),buf[4],lastLineA, False)
                    currentRefSeq.listOfAllignedSeqs.append(new_seq)
                    allignments.add(new_seq)

                lastLineA = False
                    
        f.close()
        print("done")

    #check if any block overlaps a gene, mark block if it does
    print("checking genes")
    for gene in listOfGenes:
        allignments.checkGene(gene)
    print("done")
            
            
    #mark overlapping blocks from reference genome
    #mark all of the sequences that were alligned with overlapping
    #reference blocks. Remove the marked sequences
    blocknum = 1
    for species in allignments.listOfSpecies:
        if species.name == refSpecies:
            refSpeciesObj = species
    print("looking at ref species...")
    for chromosome in refSpeciesObj.listOfChromosomes:
        for refSeq in chromosome.listOfMultiZ:
            if refSeq.Overlap == True:
                for seq in refSeq.listOfAllignedSeqs:
                    seq.Overlap = True
                chromosome.listOfMultiZ.remove(refSeq)
                
            else:
                refSeq.blockNum = blocknum
                for seq in refSeq.listOfAllignedSeqs:
                    seq.blockNum = blocknum
                blocknum +=1
    print("done")
    #check all the species for blocks over lapping with genes
    #if they overlap with genes just remove that one block
    #even if they are a reference sequence.
    print("looking at geneOverlap...")
    for species in allignments.listOfSpecies:
        for chromosome in species.listOfChromosomes:
            for seq in chromosome.listOfMultiZ:
                if seq.geneOverlap == True:
                    chromosome.listOfMultiZ.remove(seq)


    print("done")
    #go through all non-reference sequences and remove those that were associated with
    #an overlapping reference sequence.
    print("removing sequence overlaps...")
    for species in allignments.listOfSpecies:
        if species.name != refSpecies:
            for chromosome in species.listOfChromosomes:
                for seq in chromosome.listOfMultiZ:
                    if seq.Overlap == True:
                        chromosome.listOfMultiZ.remove(seq)
    print("done")
    print("done")
    return allignments

def writeBED(filePath, allignments):

    #check if there is a bed folder, if there is delete it
    if 'bed' in listdir(filePath):
        subprocess.call('rm -r '+filePath+'/bed', shell=True)
    #make a new directory /bed to store the bed files
    subprocess.call('mkdir '+filePath+'/bed', shell=True)
        
    i = 0
    for species in allignments.listOfSpecies:
        try:
            fileName = filePath+'/bed/'+species.name+'.bed'
            f = open(fileName, 'w')
        except IOError:
            raise Exception("Cannot open {}".format(fileName))
        else:
            f.seek(0)
            for chromo in species.listOfChromosomes:
                #print("Num MultiZ in {}: {}".format(chromo.name,len(chromo.listOfMultiZ)))
                previous = None
                for i in range(len(chromo.listOfMultiZ)):#The multiZs are sorted by block num as they are added
                    multiZ = chromo.listOfMultiZ[i]
                    if i -1 >= 0:
                        multiZ.fivePrime = chromo.listOfMultiZ[i-1].blockNum
                    if i +1 <len(chromo.listOfMultiZ):
                        multiZ.threePrime = chromo.listOfMultiZ[i+1].blockNum
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
                            multiZ.chromosome, multiZ.species+"_"+str(multiZ.blockNum),\
                            multiZ.s, multiZ.getEndPos(),multiZ.strand, multiZ.fivePrime,\
                            multiZ.threePrime))
                            #creates a tab separated line in detailing the block

            f.close()
            i += 1

    return i #return number of files written to

def writeGenes(filePath, allignments, genes):

    fileName = filePath+'/genes/genes.bed'

    #check if there is a bed folder, if there is delete it
    if 'genes' in listdir(filePath):
        subprocess.call('rm -r '+filePath+'/genes', shell=True)
    #make a new directory /genes to store the gene files in bed format
    subprocess.call('mkdir '+filePath+'/genes', shell=True)

    try:
        f = open(fileName, 'w')
    except IOError:
        raise Exception("Cannot open {}".format(fileName))
    else:
        f.seek(0)
        i = 1
        for gene in genes:
            gene.fivePrime, gene.threePrime = allignments.getAdjBlocks(gene) or (None, None)
            #print(allignments.getAdjBlocks(gene))
            #value = allignments.getAdjBlocks(gene)
            #print("value: {}".format(value))
            #gene.fivePrime = value
            print("gene value: {}".format(gene.fivePrime))
            
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
                    gene.chromosome, i, gene.s, gene.getEndPos(),\
                    gene.strand, gene.fivePrime, gene.threePrime))

            i+=1
        f.close()
        #sort file by chromosome them start of gene
        subprocess.call("sort -k1,1 -k2,2n "+fileName+" > "+filePath+"/genes/genes_sort.bed", shell=True)
    return


                        
def maf2bed(fileNames, outputFile, listOfGenes):
    allignments = parseMAF(fileNames, listOfGenes)
    #bedFiles = [(f[:-3]+'bed') for f in fileNames]
    writeBED(outputFile, allignments)
    writeGenes(outputFile, allignments, listOfGenes)


#test.parse("/scr/k61san/trnaevo/MultiZ/19mammals_withHuman/chrUn_GL000195v1.maf.gz")

#if this python file is called from the command line
if __name__ == "__main__":
    geneList = parseInfernal('../test/Insects/infernalIn')
    multiSeqFiles = [join('../test/Insects/mafUnzipped',f) for f in listdir('../test/Insects/mafUnzipped') if isfile(join('../test/Insects/mafUnzipped',f))]

    allignments = parseMAF(multiSeqFiles, geneList)

    #print("writing block bed files...")
    #writeBED("../test", allignments)
    #print("done")
    print("writing gene bed files...")
    writeGenes("../test", allignments, geneList)
    print("done")
