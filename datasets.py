# Imports
from collections import defaultdict
from termcolor import colored
from json import dumps
import sys
import re
import os



# Color macros and formating!
tabulation = f'\t'
taxonLeftFormat = 40
statusCenterFormat = 20
red = lambda x: colored(x, 'red')
green = lambda x: colored(x, 'green')
blue = lambda x: colored(x, 'blue')
yellow = lambda x: colored(x, 'yellow')
magenta = lambda x: colored(x, 'magenta')
cyan = lambda x: colored(x, 'cyan')



def printCollection(taxon, numberOfOrganisms, status):
    if status:
        print(f'{tabulation}{yellow(taxon):<40} {green("COLLECTED"):^20} | {blue(numberOfOrganisms)} Organism{"s" if numberOfOrganisms > 1 else ""}')
    else:
        print(f'{tabulation}{yellow(taxon):<40} {red("FAILED"):^20} | {blue(numberOfOrganisms)} Organisms')

def sizeOf(num, suffix="B"):
    for unit in ("", "K", "M", "G", "T", "P", "E", "Z"):
        if abs(num) < 1024.0:
            return [float(f'{num:3.1f}'), f'{unit}{suffix}']
        num /= 1024.0
    return [float(f'{num:3.1f}'), f'Y{suffix}']

def checkSize(species, speciesName, bigfile, verbose=1):
    calculated, calculatedUnit = species[speciesName]['filesize-unit']
    
    if not bigfile:
        realText = os.popen('du -sh "out.txt"').read()
        try:
            realMatch = re.search('([0-9]{0,3}[,]{0,1}[0-9]*)([KMGTPEY])', realText)
            real, realUnit = realMatch.group(1, 2)
        except:
            print(tabulation + red('COULD NOT FIND FILE!'))
            real, realUnit = 0, ''
    else:
        print(tabulation + red('Big file detected! Only showing calculated size!'))
        real, realUnit = 0, ''

    if verbose:
        print(f'{tabulation}{blue(calculated).replace(".", ","):<17}{magenta(calculatedUnit)} \tCalculated')
        if not bigfile:
            print(f'{tabulation}{blue(real):<17}{magenta(realUnit + 'B')} \tReal')

    return real, calculated

def pathCheckCreation(path, returnPath=False):
    if not os.path.isdir(path):
        os.mkdir(path)
    os.chdir(path)

    if returnPath:
        return os.getcwd()
    return



# Creation of the Genomes path
path = os.getcwd()
path = pathCheckCreation(path + '/' + sys.argv[1], returnPath=True)
# print(path)



# Printing dictionaries
def pretty(value, sort_keys=True, indent=4):
    print(dumps(value, sort_keys=sort_keys, indent=indent))




# Collection of organism data (accession and assembly level)
def collectInfo(taxons, verbose=1, debug=0):
    species = defaultdict(lambda: {'accession': None, 'level': None})
    for taxon in taxons:
        command = f'datasets summary genome taxon "{taxon}" --reference --assembly-source RefSeq'
        shell = os.popen(command)
        summaryRead = shell.read()
        summary = eval(summaryRead.replace('true', '"true"'))
        shell.close()

        try:
            for report in summary['reports']:
                speciesName = report['organism']['organism_name']
                species[speciesName]['accession'] = report['accession']
                species[speciesName]['level'] = report['assembly_info']['assembly_level']
                species[speciesName]['filesize-unit'] = sizeOf(int(report["assembly_stats"]["total_sequence_length"]))
            
            if verbose:
                numberOfOrganisms = len(summary['reports'])
                printCollection(taxon, numberOfOrganisms, 1)
        except:
            printCollection(taxon, 0, 0)

    if verbose:
        print()

    if debug:
        print('Species = ', end='')
        pretty(species)
        print()

    return species

# Downloading genomes
def downloadGenomes(organisms, verbose=1, checkSizes=0, progressbar=0, sizeLimit=15):
    global path

    for organismName in organisms:
        organism = organisms[organismName]
        accession = organism['accession']
        size, unit = organism['filesize-unit']
        arguments = ''

        if (checkSizes or progressbar or verbose):
            print(f'{tabulation}{yellow(organismName + f" ({accession})")}:')

        bigfile = (not (unit in ['B', 'KB', 'MB']) or ((size > sizeLimit) and not (unit in ['B', 'KB'])))
        if bigfile:
            print(f'{tabulation}{red("The estimated filesize is bigger than 10 MB!")}')
            print(f'{tabulation}{cyan("DOWNLOADING DEHYDRATED FILES")}')
            arguments += ' --dehydrated'

        if not progressbar:
            arguments += ' --no-progressbar'

        pathCheckCreation(path + '/' + accession)
        os.system(f'touch "{organismName}.txt"')
        os.system(f'datasets download genome accession {accession} --filename "out.zip"{arguments}')
        
        if not bigfile:
            os.system(f'unzip -p "out.zip" "ncbi_dataset/data/G*" > "out.txt"')
        
        if verbose:
            print(f'{tabulation}{green("DOWNLOADED")}')

        if checkSizes:
            sizes = checkSize(organisms, organismName, bigfile)

        if verbose:
            print()
        
        os.chdir(path)
            
    if verbose:
        print()

taxons = [
    # 'Pseudomonas fluorescens',
    'Mycoplasmoides pneumoniae',
    'Escherichia coli',
    # 'Orientia tsutsugamushi',
    # 'Clostridium botulinum',
    # 'Clostridium tetani',
    # 'Clostridium difficile',
    # 'Yersinia pestis',
    # 'Chromatium Okenii',

    'Vaccinium macrocarpon',
    'Orchidaceae',

    # 'canis lupus familiaris'
]
taxons.sort()

species = collectInfo(taxons)
downloadGenomes(species, checkSizes=1, progressbar=1)