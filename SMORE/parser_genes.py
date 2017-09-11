import subprocess
#import argparse
from os import listdir
from os.path import isfile, join
from dataStructures_lessMem import Gene

'''
files given to the gene parser should be:
  -in the modified BED format below.
  -end with .bed (or .bed.gz, .bed.Z, .bed.bz2)
  -species in list must match those in the MAF files given

File Format:
  -Heading (optional):
    -line must start with a '#'
  -Body:
    -one line = one gene
    -tab separated elements
    -chromosome     species     startCoord     endCoord     +or-     sequence     secondary_structure(optional)
'''
'''
parser = argparse.ArgumentParser()
#parser.add_argument("input_file", type=str)
parser.add_argument("outputDir", type=str)

args = parser.parse_args()
'''

def parseGenes(input_file, speciesList):
    listOfGenes = list()
    geneNumber = 0
    emptySpecList = list()
    
    f = catchExceptions(input_file)
    #input file sorted by species
    
    for line in f:
        if line[0] in '# \n':
            pass
        else:
            geneNumber += 1
            lineList = line.split('\t')
            if lineList[3] not in emptySpecList:
                emptySpecList.append(lineList[3])
                #print("species in specieslist: {}".format(lineList[3]))
                listOfGenes.append(list())
                if lineList[3] not in speciesList:
                    speciesList.append(lineList[3])
            if len(lineList) != 10:
                raise Exception("Input file does not have all of the fields. File shoud include 10 tab\nseparated elements. All elements not present should have an 'NA' instead. Wrong line: '{}'".format(line))
            #                                              panTro       chr3         84545             85
            listOfGenes[len(listOfGenes)-1].append(Gene(lineList[3], lineList[0], int(lineList[1]), getGeneLength(int(lineList[1]), int(lineList[2])),\
                                                        lineList[4], geneNumber, lineList[8], lineList[7], ownGene=True, _type=lineList[5], pseudoGene=lineList[6], comment=lineList[9]))
            #                                              +            32          ((___>><<) ATTCGTAGCAT               tRNA                False                   NA

            
    f.close()
    if f.name != input_file:
        subprocess.call('rm '+f.name, shell=True)

    if len(listOfGenes) == 0:
        raise Exception("no genes found!")

    return listOfGenes, speciesList

def catchExceptions(fileName):
    '''
    checks if file is in correct format( '.bed', '.bed.gz', '.bed.Z', '.bed.bz2')
    checks if file is zipped
    checks if file can be opened
    '''
    
    if not isfile(fileName):
        raise Exception("'{}' is not a file.".format(fileName))
    
    if fileName.endswith('.bed'):
        pass
    elif fileName.endswith(('.bed.gz','.bed.Z','.bed.bz2')):
        unzippedname = fileName.rstrip('.gz.Z.bz2')
        subprocess.call('zcat '+fileName+' > '+unzippedname)
        fileName = unzippedname
    else:
        raise Exception("'{}' not of bed format. Use a .bed file or gziped .bed file('.bed.gz','.bed.Z','.bed.bz2')".format(fileName))

    try:
        f = open(fileName, 'r')
    except IOError:
        raise Exception("cannot open {}".format(fileName))

    return f

def getGeneLength(start, end):
    '''
    returns the length of the gene
    
    if we are on the + strand we read left to right so end is bigger so subtract start from end.
    If we are on the - strand we read right to left so start is bigger so subtract end from start.
    '''
    if end > start:
        return end - start
    return start - end


if __name__ == '__main__':

    geneObjects = parseGenes('../textFiles/test/','../textFiles/outPut')
    print("len: {}".format(len(geneObjects)))
    print(geneObjects)
