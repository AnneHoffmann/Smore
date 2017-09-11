
import subprocess
from os import listdir
from os.path import isfile, join
from parser_infernal import parseInfernal

def endSlash(args):
    '''
    takes a list of Args. 
    returns a tuple of args with '/' at the end
    '''

    returnlist = list()
    for arg in args:
        if not arg.endswith('/'):
            returnlist.append(arg+'/')
        else:
            returnlist.append(arg)
    if len(returnlist) == 1:
        return returnlist[0]
    return returnlist

def makeFullPath(pathToRepo):
    '''
    Takes the path to main.pl
    Convertes it into the path for maf2TempBed.py
    '''
    if not pathToRepo.endswith('/pipeline/'):
        if pathToRepo.endswith('/scripts_cam/'):
            return pathToRepo
        else:
            return pathToRepo + 'pipeline/scripts_cam/'
    else:
        return pathToRepo + 'scripts_cam/'

def makeSpeciesList(_list, genomes):
    '''
    Takes an initialized species list [referenceSpecies]
    returns a full list of all the species in the genome directory
    '''
    reference = _list[0]
    for f in listdir(genomes):
        if isfile(join(genomes,f)):
            species = f.split('.')[0]
            if species != reference:
                _list.append(species)
    return _list
                
def makeOutPutDir(outPutDir):
    '''
    Determines if outPutDir exists
    if it does not: make the directory
    '''
    try:
        listdir(outPutDir)
    except FileNotFoundError:
        subprocess.call('mkdir '+outPutDir, shell=True)

def infernal(outPutDir, genomes, rnaModel, incE, incT, path):
    '''
    
    '''
    if len(genomes) == 0:
        raise Exception("No genomes given. make sure all paths given are paths to directories not files")
    if 'infernalIn' in listdir(outPutDir):
        subprocess.call('rm -r '+outPutDir+'infernalIn',shell=True)
    subprocess.call('mkdir '+outPutDir+'infernalIn',shell=True)

    i = 1 #counter for number of genes searched

    if incT == None:
        infernalOptArgs = ''
    else:
        infernalOptArgs = ' --incT '+str(incT)
    if incE == None:
        pass
    else:
        infernalOptArgs = infernalOptArgs + ' --incE '+str(incE)

    for gFile in genomes:
        species = gFile.split('/')[-1].split('.')[0]
        subprocess.call(path+'/cmsearch -o '+outPutDir+'infernalIn/'+species+'.out '+infernalOptArgs+' '+rnaModel+' '+gFile, shell=True)

    geneObjs, versionInfo = parseInfernal(outPutDir+'/infernalIn')
    return geneObjs, versionInfo

if __name__ == '__main__':
    print(makeFullPath('homes/Desktop/Git/'))
    
    a = endSlash(['/homes/biertruck/cameron/Desktop/Project_May_June_2016/GitRepo/'])
    print(makeFullPath(a))
