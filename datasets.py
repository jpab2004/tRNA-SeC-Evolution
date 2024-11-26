# Imports
from collections import defaultdict, deque
import matplotlib.pyplot as plt
from termcolor import colored
from textwrap import wrap
import os, sys, glob, ast
from json import dumps
import random, pickle
import pandas as pd
import math



# import warnings
# warnings.filterwarnings( "ignore", module = ".*matplotlib.*" )



# Color macros and formating!
red = lambda x: colored(x, 'red')
green = lambda x: colored(x, 'green')
blue = lambda x: colored(x, 'blue')
yellow = lambda x: colored(x, 'yellow')
magenta = lambda x: colored(x, 'magenta')
cyan = lambda x: colored(x, 'cyan')
numberP = yellow

tabulation = f'\t\t\t\t\t'
separator = lambda: print(magenta(f'\n{tabulation}' + ''.join(['=' for _ in range(110)]) + '\n'))

pprint = lambda x: print(f'{tabulation}{x}')

ps = lambda x: print(f'{tabulation}{x} | ', end='', flush=True)
pm = lambda x: print(f'{x} | ', end='', flush=True)
pe = lambda x: print(x)

psT = lambda x: print(f'{tabulation}{x}: ', end='', flush=True)
peT = lambda x: print(x)

psTaxa = lambda: print(tabulation, end='', flush=True)
pmTaxa = lambda x: print(f'{x} -> ', end='', flush=True)

# randomColor = lambda: '#' + ''.join([random.choice('0123456789ABCDEF') for _ in range(6)])
randomColor = lambda: '#' + ''.join([random.choice('789ABCDEF') for _ in range(6)])
randomShape = lambda: random.choice(allAvailableShapes)



# Microreact colors
HH, hh = 'FF', 'BB'

colorNone = f'#{hh}{hh}{hh}'
colorRed = f'#{HH}{hh}{hh}'
colorGreen = f'#{hh}{HH}{hh}'
colorBlue = f'#{hh}{hh}{HH}'

globalTaxonColors = defaultdict(lambda: defaultdict(lambda: None))
globalTaxonColors['None'] = colorNone

globalTaxonColors['superkingdom']['Archaea'] = colorGreen
globalTaxonColors['superkingdom']['Bacteria'] = colorRed
globalTaxonColors['superkingdom']['Eukaryota'] = colorBlue

globalTaxonColors['kingdom']['Bacillati'] = f'#F200B2'
globalTaxonColors['kingdom']['Metazoa'] = f'#00DEF2'
globalTaxonColors['kingdom']['Methanobacteriati'] = f'#DEF200'
globalTaxonColors['kingdom']['Promethearchaeati'] = f'#70218A'
globalTaxonColors['kingdom']['Viridiplantae'] = f'#3A943A'



# Micoreact shapes
allAvailableShapes = [
    "circle", "diamond", "dot", "heptagon", "heptagon-inverted", "heptagram", "heptagram-inverted", "hexagon", "hexagram", "octagon", "octagram", "pentagon",
    "pentagon-inverted", "pentagram", "pentagram-inverted", "plus", "square", "star", "tetragram", "triangle", "triangle-inverted", "triangle-right",
    "triangle-left", "chevron", "double-chevron", "chevron-inverted", "double-chevron-inverted", "chevron-right", "double-chevron-right", "chevron-left",
    "double-chevron-left", "wye", "wye-inverted"
]
shapeNone = 'cross'

globalTaxonShapes = defaultdict(lambda: defaultdict(lambda: None))
globalTaxonShapes['None'] = shapeNone

globalTaxonShapes['superkingdom']['Archaea'] = 'square'
globalTaxonShapes['superkingdom']['Bacteria'] = 'circle'
globalTaxonShapes['superkingdom']['Eukaryota'] = 'diamond'



# Microreact infos to add
globalMetadataInfos = {
    'colored': ['superkingdom', 'kingdom', 'phylum'],
    'shaped': []
}



# Consts
globalGenomesPath = None
createdGenomesPath = False

globalSpeciesFile = 'species'

globalFetchFile = None
createdFetchFile = False

globalReadyFile = None
createdReadyFile = False

globalDetectedFile = None
createdDetectedFile = False

globalTaxonomyFile = None
createdTaxonomyFile = False

globalCollectedRSSUFile = None
createdCollectedRSSUFile = False

globalProcessedFile = None
createdProcessedFile = False
globalProcessedRSSUFile = None
createdProcessedRSSUFile = False

globalMetadataFile = None
createdMetadataFile = False

globalAlignFile = None
createdAlignFile = False
globalAlignRSSUFile = None
createdAlignRSSUFile = False

globalTRNACount = None

globalTaxonLevels = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
globalTaxonLevelsSequence = {
    'superkingdom': ['superkingdom'],
    'kingdom': ['superkingdom', 'kingdom'],
    'phylum': ['superkingdom', 'kingdom', 'phylum'],
    'class': ['superkingdom', 'kingdom', 'phylum', 'class'],
    'order': ['superkingdom', 'kingdom', 'phylum', 'class', 'order'],
    'family': ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family'],
    'genus': ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus'],
    'species': ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
}

globalRSSUFileHeader = [
    'primaryAccession', 'insdcAccession', 'location', 'cutoffHead', 'cutoffTail', 'type', 'evidence', 'evidenceDescription',
    'sequenceQuality', 'alignmentQuality', 'subset', 'ncbiTaxId', 'reference', 'classification', 'silvaUri', 'sequence'
]

globalTaxonTranslation = {
    'superkingdom': 'Domínio',
    'kingdom': 'Reino',
    'phylum': 'Filo',
    'class': 'Classe',
    'order': 'Ordem',
    'family': 'Família',
    'genus': 'Gênero',
    'species': 'Espécie'
}



def printCollection(taxon, numberOfOrganisms, status, popularName=None, summaryRead=None):
    pop = f'({popularName})' if popularName != None else ''

    if status == 0:
        print(f'{tabulation}{yellow(taxon):<40} {cyan(pop):<30} {red("Failed - No reference genome"):>40} | {numberP(numberOfOrganisms)} Organisms')
    elif status == 1:
        print(f'{tabulation}{yellow(taxon):<40} {cyan(pop):<30} {green("Collected"):>40} | {numberP(numberOfOrganisms)} Organism{"s" if numberOfOrganisms > 1 else ""}')
    elif status == 2:
        print(f'{tabulation}{yellow(taxon):<40} {cyan(pop):<30} {magenta("No sequenced genome data"):>40} |')
    elif status == 3:
        print(f'{tabulation}{red("NEW DATASETS CLIENT VERSION!!!"):^120}')
        sys.exit()
    elif status == 4:
        print(f'{tabulation}{yellow(taxon):<40} {cyan(pop):<30} {red("Collection error"):>40} | {red(summaryRead.replace("\n", " "))}')
    else:
        print(f'{tabulation}{red("ERROR! INVALID COLLECTION STATUS!")}')


def sizeOf(num, suffix="B"):
    for unit in ("", "K", "M", "G", "T", "P", "E", "Z"):
        if abs(num) < 1024.0:
            return [float(f'{num:3.1f}'), f'{unit}{suffix}']
        num /= 1024.0
    return [float(f'{num:3.1f}'), f'Y{suffix}']

def checkPath(path, returnPath=False):
    if not os.path.isdir(path):
        os.mkdir(path)

    if returnPath:
        return path

def pathCheckCreation(path, returnPath=False, create=True):
    if create:
        if not os.path.isdir(path):
            os.mkdir(path)
    
    os.chdir(path)

    if returnPath:
        return os.getcwd()
    return

def checkTaxonomyRSSU(organismName, allInfos, taxonomyInfos):
    global globalTaxonLevels

    for globalTaxon in globalTaxonLevels[::-1]:
        if not globalTaxon in taxonomyInfos[organismName]['taxonomy']: continue
        specificTaxon = taxonomyInfos[organismName]['taxonomy'][globalTaxon]['name']

        for entry in allInfos:
            entryInfos = allInfos[entry]
            entryClassification = entryInfos['classification']
            
            if specificTaxon in entryClassification:
                rType = entryInfos['type']
                sequence = entryInfos['sequence']
                header = f'>{organismName};{rType}'

                return True, header, sequence
    
    return False, None, None

def addColors(c1, c2):
    c1 = wrap(c1.replace('#', ''), 2)
    c2 = wrap(c2.replace('#', ''), 2)
    
    c1 = [int(c, 16) for c in c1]
    c2 = [int(c, 16) for c in c2]

    color = [int((c1[i] + c2[i])/2) for i in range(len(c1))]
    color = '#' + ''.join([f'{c:0{2}x}' for c in color])

    return color

def plotTaxonAnalysis(taxonAnalysis, levels, sequential=False, totalColor='#999999', colorWeight='#222222', unique=False, show=False):
    def annotateAx(taxon, found, total, percentage, i, c1, c2, hatch):
        edgeColor, lineWidth = '#000000', 1

        if len(taxons) / fig.get_figwidth() < 1.5:
            ax.bar(taxon, found, color=c2, bottom=0, edgecolor=edgeColor, linewidth=lineWidth, hatch=hatch)
            ax.bar(taxon, total-found, color=c1, bottom=found, edgecolor=edgeColor, linewidth=lineWidth)
            
            totalY = 1.2*found
            ax.annotate(f'{percentage}%', (i, totalY), ha='center', va='center', rotation=0)
            return False
        else:
            ax.barh(taxon, found, color=c2, left=0, edgecolor=edgeColor, linewidth=lineWidth, hatch=hatch)
            ax.barh(taxon, total-found, color=c1, left=found, edgecolor=edgeColor, linewidth=lineWidth)

            totalY = total + .1*total
            if percentage > 95:
                ax.annotate('*', (totalY, i), ha='center', va='center', rotation=0)
                return True
            if percentage <= 0:
                ax.annotate('*', (totalY, i), ha='center', va='center', rotation=0, color='#ff0000')
                return True

    global globalTaxonColors
    figs, axes = [], []

    for level in levels:
        annotated = False

        taxons = list(taxon for taxon in taxonAnalysis[level] if not taxon in ['found', 'total', 'percentage', 'counted'])
        taxons = sorted(taxons, key=lambda x: taxonAnalysis[level][x]['total'] + 1e-3*taxonAnalysis[level][x]['found'], reverse=True)
        hatch = f'/'*math.ceil((len(taxons)+1)/20)
        
        if len(taxons) / 12 < 1.5:
            fig, ax = plt.subplots(1, 1, figsize=(12, 6), dpi=150)
        else:
            fig, ax = plt.subplots(1, 1, figsize=(12, 20), dpi=150)

        found = taxonAnalysis[level]['found']
        total = taxonAnalysis[level]['total']
        percentage = taxonAnalysis[level]['percentage']
        annotated = True if annotateAx('Total', found, total, percentage, 0, totalColor, addColors(totalColor, colorWeight), hatch) else annotated
        

        for i, taxon in enumerate(taxons, 1):
            if not taxon in globalTaxonColors[level]:
                globalTaxonColors[level][taxon] = randomColor()
            color = globalTaxonColors[level][taxon]
            colorFound = addColors(color, colorWeight)

            found = taxonAnalysis[level][taxon]['found']
            total = taxonAnalysis[level][taxon]['total']
            percentage = taxonAnalysis[level][taxon]['percentage']
            annotated = True if annotateAx(taxon, found, total, percentage, i, color, colorFound, hatch) else annotated

        xLabel = globalTaxonTranslation[level]
        orgRaw = 'species' if unique else 'organisms'
        org = 'espécies' if unique else 'organismos'

        if len(taxons) / fig.get_figwidth() < 1.5:
            ax.set_title(f'Número de tRNA-SeC por taxon ({org})', fontweight='bold', fontsize=25)
            ax.set_yscale('log')
        
            ax.set_xlabel(xLabel, fontweight='bold', fontsize=20)
            ax.set_ylabel(f'Número de tRNA-SeC (log)', fontweight='bold', fontsize=15)
            ax.grid(ls='--', alpha=.7, axis='y')
            ax.bar(0, 0, color='#ffffff', hatch='///', edgecolor='#000000', label='Contêm tRNA-SeC')

            ax.set_xticks(ax.get_xticks(), labels=ax.get_xticklabels(), fontsize=20, fontweight='bold')
            ax.set_yticks(ax.get_yticks(), labels=ax.get_yticklabels(), fontsize=15, fontweight='bold')

            ax.set_ylim([.5, .1*ax.get_ylim()[1]])
            ax.set_xlim([-1, len(taxons)+1])
        else:
            ax.set_title(f'Número de tRNA-SeC por taxon ({org})', fontweight='bold', fontsize=25, x=.25)
            ax.set_xscale('log')

            ax.set_ylabel(xLabel, fontweight='bold', fontsize=30)
            ax.set_xlabel(f'Número de tRNA-SeC (log)', fontweight='bold', fontsize=25)
            ax.grid(ls='--', alpha=.7, axis='x')
            ax.barh(0, 0, color='#ffffff', hatch='///', edgecolor='#000000', label='Contêm tRNA-SeC')

            ax.set_xticks(ax.get_xticks(), labels=ax.get_xticklabels(), fontsize=20, fontweight='bold')
            ax.set_yticks(ax.get_yticks(), labels=ax.get_yticklabels(), fontsize=15, fontweight='bold')

            ax.set_xlim([.5, .1*ax.get_xlim()[1]])
            ax.set_ylim([-1, len(taxons)+1])

        ax.legend()
        if annotated:
            ax.plot([], [], label='*: 0% do taxon', color='#ffffff')
            ax.plot([], [], label='*: > 95% do taxon', color='#ffffff')

            handles, labels = ax.get_legend_handles_labels()
            order = [2,0,1]
            legTxts = ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=20).get_texts()
            legTxts[1].set_color('#ff0000')

        fig.tight_layout()
        if ((show) and (sequential)):
            plt.show()

        pathCheckCreation(f'{globalGenomesPath}/figs')
        fig.savefig(f'{globalGenomesPath}/figs/{level}_{orgRaw}.png', format='png')

        figs.append(fig)
        axes.append(ax)
    
    if ((show) and (not sequential)):
        plt.show()

    os.system(f'cd {globalGenomesPath}')
    
    return



# Creation of files
def createGenomesPath(__genomesPath='Genomes/', suppress=False):
    global globalGenomesPath, createdGenomesPath

    __path = os.getcwd()
    __path = pathCheckCreation(__path + '/' + __genomesPath, returnPath=True, create=not suppress)
    createdGenomesPath = True

    globalGenomesPath = __path

def createFetchFile(__fetchFile='fetch', suppress=False):
    global globalFetchFile, createdFetchFile
    
    if not suppress: os.system(f'> {__fetchFile}.data')
    createdFetchFile = True

    globalFetchFile = os.popen(f'readlink -f {__fetchFile}.data').read()[:-1]

def createReadyFile(__readyFile='ready', suppress=False):
    global globalReadyFile, createdReadyFile
    
    if not suppress: os.system(f'> {__readyFile}.data')
    createdReadyFile = True

    globalReadyFile = os.popen(f'readlink -f {__readyFile}.data').read()[:-1]

def createDetectedFile(__detectedFile='detected', suppress=False):
    global globalDetectedFile, createdDetectedFile
    
    if not suppress: os.system(f'> {__detectedFile}.data')
    createdDetectedFile = True

    globalDetectedFile = os.popen(f'readlink -f {__detectedFile}.data').read()[:-1]

def createTaxonomyFile(__taxonomyFile='taxonomy', suppress=False):
    global globalTaxonomyFile, createdTaxonomyFile
    
    if not suppress: os.system(f'> {__taxonomyFile}.data')
    createdTaxonomyFile = True

    globalTaxonomyFile = os.popen(f'readlink -f {__taxonomyFile}.data').read()[:-1]

def createCollectedRSSUFile(__rSSUFile='collectedRSSU', suppress=False):
    global globalCollectedRSSUFile, createdCollectedRSSUFile
    
    if not suppress: os.system(f'> {__rSSUFile}.data')
    createdCollectedRSSUFile = True

    globalCollectedRSSUFile = os.popen(f'readlink -f {__rSSUFile}.data').read()[:-1]

def createProcessedFile(__processedFile='processed', suppress=False):
    global globalProcessedFile, createdProcessedFile
    
    if not suppress: os.system(f'> {__processedFile}.fna')
    createdProcessedFile = True

    globalProcessedFile = os.popen(f'readlink -f {__processedFile}.fna').read()[:-1]

def createProcessedRSSUFile(__processedRSSUFile='processedRSSU', suppress=False):
    global globalProcessedRSSUFile, createdProcessedRSSUFile
    
    if not suppress: os.system(f'> {__processedRSSUFile}.fna')
    createdProcessedRSSUFile = True

    globalProcessedRSSUFile = os.popen(f'readlink -f {__processedRSSUFile}.fna').read()[:-1]

def createAlignFile(__alignFile='align', suppress=False):
    global globalAlignFile, createdAlignFile
    
    if not suppress: os.system(f'> {__alignFile}.fasta')
    createdAlignFile = True

    globalAlignFile = os.popen(f'readlink -f {__alignFile}.fasta').read()[:-1]

def createAlignRSSUFile(__alignRSSUFile='alignRSSU', suppress=False):
    global globalAlignRSSUFile, createdAlignRSSUFile
    
    if not suppress: os.system(f'> {__alignRSSUFile}.fasta')
    createdAlignRSSUFile = True

    globalAlignRSSUFile = os.popen(f'readlink -f {__alignRSSUFile}.fasta').read()[:-1]

def createMetadataFile(__metadataFile='metadata', suppress=False):
    global globalMetadataFile, createdMetadataFile
    
    if not suppress: os.system(f'> {__metadataFile}.csv')
    createdMetadataFile = True

    globalMetadataFile = os.popen(f'readlink -f {__metadataFile}.csv').read()[:-1]

def initiate(
        __genomesPath='Genomes/',
        __fetchFile='fetch',
        __readyFile='ready',
        __detectedFile='detected',
        __taxonomyFile='taxonomy',
        __rSSUFile='rSSU',
        __processedFile='processed',
        __processedRSSUFile='processedRSSU',
        __alignFile='align',
        __alignRSSUFile='alignRSSU',
        __metadataFile='metadata',
        
        suppressDownload=False,
        suppressFetch=False,
        suppressDetected=False,
        suppressTaxonomy=False,
        suppressRSSU=False,
        suppressProcess=False,
        suppressAlign=False
    ):
    createGenomesPath(__genomesPath, suppress=suppressDownload)
    createFetchFile(__fetchFile, suppress=suppressFetch)
    createReadyFile(__readyFile, suppress=(suppressDownload and suppressFetch))
    createDetectedFile(__detectedFile, suppress=suppressDetected)
    createTaxonomyFile(__taxonomyFile, suppress=suppressTaxonomy)
    createCollectedRSSUFile(__rSSUFile, suppress=suppressRSSU)
    createProcessedFile(__processedFile, suppress=suppressProcess)
    createProcessedRSSUFile(__processedRSSUFile, suppress=suppressProcess)
    createAlignFile(__alignFile, suppress=suppressAlign)
    createAlignRSSUFile(__alignRSSUFile, suppress=suppressAlign)
    createMetadataFile(__metadataFile, suppress=suppressProcess)

    return



def addToFetchFile(fetchInfos, __fetchFile=None):
    global globalFetchFile
    
    if __fetchFile == None:
        fetchFile = globalFetchFile
    else:
        fetchFile = __fetchFile

    if not createdFetchFile:
        createFetchFile(fetchFile)
        fetchFile = globalFetchFile

    with open(fetchFile, 'a') as fileHandler:
        for name in fetchInfos:
            accession = fetchInfos[name]['accession']
            taxId = fetchInfos[name]['tax-id']
            fetchFolderPath = fetchInfos[name]['fetch-folder'].strip('\n')
            chromosomesFolderPath = fetchInfos[name]['chromosomes-folder'].strip('\n')
            popularName = fetchInfos[name]['popular-name']
            kingdom = fetchInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {taxId} > {popularName} > {fetchFolderPath} > {chromosomesFolderPath} > {kingdom}\n')

def addToReadyFile(readyInfos, __readyFile=None):
    global globalReadyFile
    
    if __readyFile == None:
        readyFile = globalReadyFile
    else:
        readyFile = __readyFile

    if not createdReadyFile:
        createReadyFile(readyFile)
        readyFile = globalReadyFile

    with open(readyFile, 'a') as fileHandler:
        for name in readyInfos:
            accession = readyInfos[name]['accession']
            taxId = readyInfos[name]['tax-id']
            chromosomesFolderPath = readyInfos[name]['chromosomes-folder'].strip('\n')
            popularName = readyInfos[name]['popular-name']
            kingdom = readyInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {taxId} > {popularName} > {chromosomesFolderPath} > {kingdom}\n')

def addToDetectedFile(detectedInfos, __detectedFile=None):
    global globalDetectedFile
    
    if __detectedFile == None:
        detectedFile = globalDetectedFile
    else:
        detectedFile = __detectedFile

    if not createdDetectedFile:
        createDetectedFile(detectedFile)
        detectedFile = globalDetectedFile

    with open(detectedFile, 'a') as fileHandler:
        for name in detectedInfos:
            accession = detectedInfos[name]['accession']
            taxId = detectedInfos[name]['tax-id']
            chromosomesFolderPath = detectedInfos[name]['chromosomes-folder'].strip('\n')
            popularName = detectedInfos[name]['popular-name']
            kingdom = detectedInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {taxId} > {popularName} > {chromosomesFolderPath} > {kingdom}\n')

def addToRSSUFile(rSSUInfos, __rSSUFile=None):
    global globalCollectedRSSUFile
    
    if __rSSUFile == None:
        rSSUFile = globalCollectedRSSUFile
    else:
        rSSUFile = __rSSUFile

    if not createdCollectedRSSUFile:
        createCollectedRSSUFile(rSSUFile)
        rSSUFile = globalCollectedRSSUFile

    with open(rSSUFile, 'a') as fileHandler:
        for name in rSSUInfos:
            accession = rSSUInfos[name]['accession']
            taxId = rSSUInfos[name]['tax-id']
            chromosomesFolderPath = rSSUInfos[name]['chromosomes-folder'].strip('\n')
            popularName = rSSUInfos[name]['popular-name']

            fileHandler.write(f'{accession} > {name} > {taxId} > {popularName} > {chromosomesFolderPath}\n')

def addToProcessedFile(processedInfos, __processedFile=None):
    global globalProcessedFile
    
    if __processedFile == None:
        processedFile = globalProcessedFile
    else:
        processedFile = __processedFile

    if not createdProcessedFile:
        createProcessedFile(processedFile)
        processedFile = globalProcessedFile

    with open(processedFile, 'a') as fileHandler:
        for name in processedInfos:
            header = processedInfos[name]['header']
            sequence = processedInfos[name]['sequence']

            fileHandler.write(f'{header}\n{sequence}\n')
        
def addToProcessedRSSUFile(processedRSSUInfos, __processedRSSUFile=None):
    global globalProcessedRSSUFile
    
    if __processedRSSUFile == None:
        processedRSSUFile = globalProcessedRSSUFile
    else:
        processedRSSUFile = __processedRSSUFile

    if not createdProcessedRSSUFile:
        createProcessedRSSUFile(processedRSSUFile)
        processedRSSUFile = globalProcessedRSSUFile

    with open(processedRSSUFile, 'a') as fileHandler:
        for name in processedRSSUInfos:
            header = processedRSSUInfos[name]['header']
            sequence = processedRSSUInfos[name]['sequence']

            fileHandler.write(f'{header}\n{sequence}\n')

def addToTaxonomyFile(taxonInfos, __taxonomyFile=None):
    global globalTaxonLevels, globalTaxonomyFile

    if __taxonomyFile == None:
        taxonomyFile = globalTaxonomyFile
    else:
        taxonomyFile = __taxonomyFile

    if not createdTaxonomyFile:
        createTaxonomyFile(taxonomyFile)
        taxonomyFile = globalTaxonomyFile

    with open(taxonomyFile, 'a') as fileHandler:
        for name in taxonInfos:
            accession = taxonInfos[name]['accession']
            taxId = taxonInfos[name]['tax-id']
            chromosomesFolderPath = taxonInfos[name]['chromosomes-folder'].strip('\n')
            popularName = taxonInfos[name]['popular-name']

            commandParts = []
            for level in taxonInfos[name]['taxonomy']:
                if 'id' in level: continue
                taxon, taxonId = taxonInfos[name]['taxonomy'][level], taxonInfos[name]['taxonomy'][f'{level}-id']
                commandParts.append(f'{level}<{taxon}<{taxonId}')

            command = ' | '.join(commandParts)
            fileHandler.write(f'{accession} > {name} > {taxId} > {popularName} > {chromosomesFolderPath} > {command}\n')

def addToMetadataFile(metadataInfos, __metadataFile=None):
    global globalTaxonLevels, globalTaxonColors, globalMetadataFile, globalMetadataInfos
    
    coloredTaxons = globalMetadataInfos['colored']
    shapedTaxons = globalMetadataInfos['shaped']

    taxonsTextGeneral = ';'.join(globalTaxonLevels)
    coloredTaxonsTextGeneral = ';'.join([f'{cTaxon}__color' for cTaxon in coloredTaxons])
    shapedTaxonsTextGeneral = ';'.join([f'{sTaxon}__shape' for sTaxon in shapedTaxons])

    if __metadataFile == None:
        metadataFile = globalMetadataFile
    else:
        metadataFile = __metadataFile

    if not createdMetadataFile:
        createMetadataFile(metadataFile)
        metadataFile = globalMetadataFile

    with open(metadataFile, 'a') as fileHandler:
        firstLine = f'ID;taxID;tRNA-SeC number;rRNA type;Mitochondrial rRNA;score;{taxonsTextGeneral};{coloredTaxonsTextGeneral};{shapedTaxonsTextGeneral}\n'
        fileHandler.write(firstLine)

        for identification in metadataInfos:
            taxId = metadataInfos[identification]['tax-id']
            taxonomy = metadataInfos[identification]['taxonomy']
            score = metadataInfos[identification]['score']
            tRNANumber = metadataInfos[identification]['tRNA-number']
            rType = metadataInfos[identification]['rRNA-type']
            mitochondrial = metadataInfos[identification]['mitochondrial']

            taxons = []
            for __level in globalTaxonLevels:
                taxons.append(taxonomy[__level] if __level in taxonomy else str(None))
            taxonsText = ';'.join(taxons)

            taxonsColors = []
            for superLevel in coloredTaxons:
                if superLevel in taxonomy:
                    if not taxonomy[superLevel] in globalTaxonColors[superLevel]:
                        globalTaxonColors[superLevel][taxonomy[superLevel]] = randomColor()
                    colorToAppend = globalTaxonColors[superLevel][taxonomy[superLevel]]
                else:
                    colorToAppend = globalTaxonColors['None']    
                taxonsColors.append(colorToAppend)
            
            taxonsShapes = []
            for superLevel in shapedTaxons:
                if superLevel in taxonomy:
                    if not taxonomy[superLevel] in globalTaxonShapes[superLevel]:
                        globalTaxonShapes[superLevel][taxonomy[superLevel]] = randomShape()
                    shapeToAppend = globalTaxonShapes[superLevel][taxonomy[superLevel]]
                else:
                    shapeToAppend = globalTaxonShapes['None']
                taxonsShapes.append(shapeToAppend)

            taxonsColorText = ';'.join(taxonsColors)
            taxonsShapeText = ';'.join(taxonsShapes)

            text = f'{identification};{taxId};{tRNANumber};{rType};{mitochondrial};{score};{taxonsText};{taxonsColorText};{taxonsShapeText}\n'

            fileHandler.write(text)



# Printing dictionaries
def pretty(value, sort_keys=True, indent=4, colorOrder=[green, blue, red, cyan, yellow, magenta], colorOffset=False):
    color = deque(colorOrder)
    color.rotate(colorOffset)
    
    dump = str(repr(dumps(value, sort_keys=sort_keys, indent=indent))).split(r'\n')
    for packet in dump:
        packet = packet.strip("{},'")
        if packet != '':
            depth = packet.count('    ')

            color.rotate(depth-1)
            pcolor = color[0]
            color.rotate(- (depth-1))
            pprint(pcolor(packet[4:]))




# Collection of organism data (accession and assembly level)
def collectInfo(taxons, verbose=True, archaea=False, save=False, read=False, __speciesFile=None):
    species = defaultdict(lambda: {'accession': None})

    if __speciesFile == None:
        global globalSpeciesFile
        speciesFile = globalSpeciesFile + '.pickle'
    else:
        speciesFile = __speciesFile + '.pickle'

    separator()
    print(f'{tabulation}{magenta("Data collection starting"):^121}\n')

    if ((read) and (not save)):
        print(f'{tabulation}{green("Read flag detected") + " | " + magenta("Trying to read file"):^130}\n')

        with open(speciesFile, 'rb') as fileHandler:
            species = pickle.load(fileHandler)

        print(f'{tabulation}{(green("Read successful (") + numberP(len(list(species.keys()))) + green(" genomes)")):^140}')
        separator()
    else:
        totalIndexing = len(list(taxons.keys()))
        indexingSize = len(str(totalIndexing))
        numberOfOrganismsTotal = 0
        for i, (taxon, (popularName, kingdom)) in enumerate(sorted(taxons.items()), 1):
            command = f'datasets summary genome taxon "{taxon}"'
            command += f' --reference --assembly-source RefSeq --as-json-lines 2>&1' if (not archaea) else ' --assembly-source RefSeq --as-json-lines 2>&1'
            shell = os.popen(command)
            summaryRead = shell.read().replace('true', "'true'")[:-1]
            name = f'{i:<{indexingSize}}/{totalIndexing}. {taxon}'
            shell.close()

            if 'but no genome data is currently available for this taxon' in summaryRead:
                printCollection(name, 0, 2, popularName=popularName)
                continue
            elif 'New version of client' in summaryRead:
                printCollection('', 0, 3)
                continue

            summary = {}
            for j, summa in enumerate(summaryRead.split('\n')):
                try:
                    summary[j] = ast.literal_eval(summa)
                except Exception as e:
                    print(j, summa, e)
                    sys.exit()
            
            numberOfOrganisms = 0
            try:
                for report in summary.values():
                    try:
                        speciesName = report['organism']['organism_name']
                        species[speciesName]['tax-id'] = report['organism']['tax_id']
                        species[speciesName]['accession'] = report['accession']
                        species[speciesName]['filesize-unit'] = sizeOf(int(report["assembly_stats"]["total_sequence_length"]))
                        species[speciesName]['popular-name'] = popularName
                        species[speciesName]['kingdom'] = kingdom
                        numberOfOrganisms += 1
                        numberOfOrganismsTotal += 1
                    except Exception as e:
                        print(summary, e)

                if verbose:
                    printCollection(name, numberOfOrganisms, 1, popularName=popularName)
            except Exception as e:
                print(report, e)
                printCollection(name, 0, 0, popularName=popularName)

        print(f'\n{tabulation}{magenta(f"Data collection ended") + " | " + numberP(numberOfOrganismsTotal) + magenta(" total organisms collected"):^140}')
        separator()

    if save:
        with open(speciesFile, 'wb') as fileHandler:
            aux = dict(species)
            pickle.dump(aux, fileHandler, protocol=pickle.HIGHEST_PROTOCOL)

    return species

# Downloading genomes
def downloadGenomes(
        organisms,
        __genomesPath=None,
        __readyFile=None,
        __fetchFile=None,
        referenceRange=None,
        rangeStep=None,
        verbose=True,
        reDownload=False,
        progressbar=False,
        sizeLimit=15,
        zipFile='genome',
        fetchFolder='fetchFolder',
        chromosomesFolder='chromosomes'
    ):

    if __genomesPath == None:
        global globalGenomesPath, createdGenomesPath
        genomesPath = globalGenomesPath
    else:
        genomesPath = __genomesPath

    lost = 0
    found = 0
    skipped = 0
    fetchInfos = {}
    readyInfos = {}

    if not createdGenomesPath:
        createGenomesPath()

    if not createdReadyFile:
        createReadyFile()

    if referenceRange == None:
        organismsNamesSorted = list(organisms.keys())
    else:
        if rangeStep == None:
            rangeStep = 1

        range1 = referenceRange * (rangeStep - 1)
        range2 = referenceRange + (referenceRange * (rangeStep - 1))
        
        if range1 > len(list(organisms.keys())):
            range1 = len(list(organisms.keys())) - referenceRange
        if range2 > len(list(organisms.keys())):
            range2 = len(list(organisms.keys()))

        organismsNamesSorted = list(organisms.keys())[range1:range2]

    print(f'{tabulation}{magenta("Genomes download starting"):^120}\n')
    totalIndexing = len(organismsNamesSorted)
    indexingSize = len(str(totalIndexing))
    for i, organismName in enumerate(organismsNamesSorted, 1):
        organism = organisms[organismName]
        accession = organism['accession']
        size, unit = organism['filesize-unit']
        popularName = organism['popular-name']
        kingdom = organism['kingdom']
        name = f'{i:<{indexingSize}}/{totalIndexing}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')
        arguments = ''

        shellRecycled = os.popen(f'if [ -e {genomesPath}/{accession}/recycled.status ]; then echo "1"; else echo "0"; fi;')
        shellDownloaded = os.popen(f'if [ -e {genomesPath}/{accession}/downloaded.status ]; then echo "1"; else echo "0"; fi;')
        shellRehydrated = os.popen(f'if [ -e {genomesPath}/{accession}/rehydrated.status ]; then echo "1"; else echo "0"; fi;')
        shellNoChromosome = os.popen(f'if [ -e {genomesPath}/{accession}/noChromosome.status ]; then echo "1"; else echo "0"; fi;')

        bigfile = (not (unit in ['B', 'KB', 'MB']) or ((size > sizeLimit) and not (unit in ['B', 'KB'])))

        if (progressbar or verbose):
            print(f'{tabulation}{yellow(name)}:')

        if reDownload:
            shellRemoveRecyled = os.popen(f'rm "{genomesPath}/{accession}/recycled.status"')
            _ = shellRemoveRecyled.read()
            shellRemoveRecyled.close()
        elif int(shellRecycled.read()):
            found += 1
            if (progressbar or verbose):
                ps(red('Organism already analysed and recycled'))
                pe(magenta('Skipping organism'))
                print()

            readyInfos[organismName] = organisms[organismName]
            readyInfos[organismName]['popular-name'] = organism['popular-name']
            readyInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath}/{accession + "/" + chromosomesFolder}').read()
            readyInfos[organismName]['kingdom'] = organism['kingdom']
            continue

        if int(shellNoChromosome.read()):
            if (progressbar or verbose):
                ps(red('Organism already marked as not having complete chromosomes'))
                pe(magenta('Skipping organism'))
                print()

            lost += 1
            continue

        if bigfile:
            if int(shellRehydrated.read()):
                found += 1
                if (progressbar or verbose):
                    ps(green('Already rehydrated'))
                    pe(magenta('Skipping organism'))
                    print()

                readyInfos[organismName] = organisms[organismName]
                readyInfos[organismName]['popular-name'] = organism['popular-name']
                readyInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath}/{accession + "/" + chromosomesFolder}').read()
                readyInfos[organismName]['kingdom'] = organism['kingdom']
                continue
            elif int(shellDownloaded.read()):
                if (progressbar or verbose):
                    ps(green('Already downloaded'))
                    pe(magenta('Skipping organism'))
                    print()

                fetchInfos[organismName] = organisms[organismName]
                fetchInfos[organismName]['fetch-folder'] = os.popen(f'readlink -f {genomesPath}/{accession + "/" + fetchFolder}').read()
                fetchInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath}/{accession + "/" + chromosomesFolder}').read()
                fetchInfos[organismName]['popular-name'] = organism['popular-name']
                continue

            ps(red(f'The estimated filesize is bigger than {sizeLimit} MB!'))
            pe(magenta('Downloading dehydrated files'))
            arguments += ' --dehydrated'
        elif int(shellDownloaded.read()):
            found += 1
            if (progressbar or verbose):
                ps(green('Already downloaded'))
                pe(magenta('Skipping organism'))
                print()

            readyInfos[organismName] = organisms[organismName]
            readyInfos[organismName]['popular-name'] = organism['popular-name']
            readyInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath}/{accession + "/" + chromosomesFolder}').read()
            readyInfos[organismName]['kingdom'] = organism['kingdom']
            continue

        if not progressbar:
            arguments += ' --no-progressbar'

        if kingdom == 'E':
            arguments += ' --chromosomes all'

        try:
            pathCheckCreation(genomesPath + '/' + accession)

            ps(green('Download folder found'))
            pm(cyan('Attempting download'))

            shellZipFolder = os.popen(f'if [ -e {genomesPath}/{accession}/{zipFile}.zip ]; then echo "1"; else echo "0"; fi;')
            if int(shellZipFolder.read()):
                if verbose:
                    pm(green('Zip folder found'))
                    pm(magenta('Extracting'))
            else:
                shellTouch = os.popen(f'> "{organismName.replace("/", "_")}.name"')
                _ = shellTouch.read()
                shellTouch.close()

                shellDownload = os.popen(f'datasets download genome accession {accession} --filename "{zipFile}.zip"{arguments} 2>&1')
                result = shellDownload.read()
                shellDownload.close()

                if 'connection reset by peer' in result:
                    if verbose:
                        pm(red('Connection reseted'))
                        pe(magenta('Skipping organism'))
                        print()
                        skipped += 1
                    continue

                if 'timeout' in result:
                    if verbose:
                        pm(red('Connection timeout'))
                        pe(magenta('Skipping organism'))
                        skipped += 1
                        print()
                    continue

                if 'Error' in result:
                    if verbose:
                        pm(red('Error'))
                        pe(magenta(result))
                        skipped += 1
                        print()
                    continue
        except KeyboardInterrupt:
            sys.exit()
        except Exception as e:
            print(tabulation + red('ERROR:') + str(e))
            sys.exit()
        

        try:
            if not bigfile:
                checkPath(genomesPath + '/' + accession + '/' + chromosomesFolder)
                folder = os.popen(f'readlink -f {genomesPath}/{accession + "/" + chromosomesFolder}').read().strip('\n')

                zipShell = os.popen(f'unzip -p "{zipFile}.zip" "ncbi_dataset/data/G*" > "{folder}/chromosome.fna"')
                _ = zipShell.read()
                zipShell.close()

                readyInfos[organismName] = organisms[organismName]
                readyInfos[organismName]['popular-name'] = organism['popular-name']
                readyInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath}/{accession + "/" + chromosomesFolder}').read()
                readyInfos[organismName]['kingdom'] = organism['kingdom']
            else:
                checkPath(genomesPath + '/' + accession + '/' + chromosomesFolder)

                zipShell = os.popen(f'unzip -o "{zipFile}.zip" -d "{fetchFolder}"')
                _ = zipShell.read()
                zipShell.close()

                fetchInfos[organismName] = organisms[organismName]
                fetchInfos[organismName]['fetch-folder'] = os.popen(f'readlink -f {genomesPath}/{accession + "/" + fetchFolder}').read()
                fetchInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath}/{accession + "/" + chromosomesFolder}').read()
                fetchInfos[organismName]['popular-name'] = organism['popular-name']
                fetchInfos[organismName]['kingdom'] = organism['kingdom']
        except KeyboardInterrupt:
            zipShell.close()
            sys.exit()
        except Exception as e:
            print(tabulation + red('ERROR:') + str(e))
        
        os.system(f'> "downloaded.status"')
        
        if verbose:
            pe(green('Downloaded'))

        if verbose:
            print()

        os.chdir(genomesPath)
        
    addToReadyFile(readyInfos, __readyFile)
    addToFetchFile(fetchInfos, __fetchFile)

    totalTemp = len(organismsNamesSorted)
    fetchTemp = len(list(fetchInfos.keys()))
    downloadedTemp = totalTemp - fetchTemp - found - lost
    print(
        f'{tabulation}{numberP(downloadedTemp) + magenta(f" genomes downloaded") + " | " + numberP(found) + magenta(" found, ")
        + numberP(skipped) + magenta(" skipped & ") + numberP(lost) + magenta(" lost & ") + numberP(fetchTemp) + magenta(" more for rehydration"):^205}'
    )
    separator()



def downloadFetch(__fetchFile=None, __readyFile=None, verbose=True, progressbar=False):
    if __fetchFile == None:
        global globalFetchFile
        fetchFile = globalFetchFile
    else:
        fetchFile = __fetchFile

    with open(fetchFile, 'r') as fileHandler:
        fetchLines = fileHandler.readlines()
    
    fetchInfos = {}
    readyInfos = {}
    for info in fetchLines:
        splitted = info.strip('\n').split(' > ')
        fetchInfos[splitted[1]] = {
            'accession': splitted[0],
            'tax-id': splitted[2],
            'popular-name': splitted[3],
            'fetch-folder': splitted[4],
            'chromosomes-folder': splitted[5],
            'kingdom': splitted[6]
        }

    lost = 0
    found = 0
    skipped = 0
    print(f'{tabulation}{magenta(f"Genomes rehydration starting"):^120}\n')
    totalIndexing = len(list(fetchInfos.keys()))
    indexingSize = len(str(totalIndexing))
    for i, organismName in enumerate(fetchInfos, 1):
        fetchFolder = fetchInfos[organismName]['fetch-folder']
        chromosomesFolder = fetchInfos[organismName]['chromosomes-folder']
        popularName = fetchInfos[organismName]['popular-name']
        accession = fetchInfos[organismName]['accession']
        kingdom = fetchInfos[organismName]['kingdom']
        name = f'{i:<{indexingSize}}/{totalIndexing}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')
        arguments = ' --directory . --max-workers 30'

        shellRehydrated = os.popen(f'if [ -e {globalGenomesPath}/{accession}/rehydrated.status ]; then echo "1"; else echo "0"; fi;')
        if int(shellRehydrated.read()):
            found += 1
            if verbose:
                print(f'{tabulation}{yellow(name)}:')
                ps(green('Already rehydrated'))
                pe(magenta('Skipping organism'))
                print()
            continue

        if (verbose):            
            print(f'{tabulation}{yellow(name)}:')
        
        try:
            os.chdir(fetchFolder)
            if verbose:
                ps(green('Fetch folder found'))
                pm(cyan('Attempting rehydration'))
        except:
            if verbose:
                ps(red('Could not find fetch folder!'))
                pe(cyan('Skipping organism'))
                print()

            skipped += 1
            continue

        if kingdom == 'E':
            arguments += ' --match chr'
        if not progressbar:
            arguments += ' --no-progressbar'

        try:
            rehydrationShell = os.popen(f'datasets rehydrate --directory .{arguments}')
            rehydrationRead = rehydrationShell.read()
            rehydrationShell.close()

            if 'Found no files for rehydration' in rehydrationRead:
                pm(red('No Complete chromosome found'))
                pe(magenta('Skipping organism'))
                print()
                
                shellTouch = os.popen('> "../noChromosome.status"')
                _ = shellTouch.read()
                shellTouch.close()

                lost += 1

                continue

            mvFilesTemp1 = fetchFolder + '/ncbi_dataset/' + os.popen(f'datasets rehydrate{arguments} --list').read().strip('\n')
            mvFilesTemp2 = mvFilesTemp1.split('\n')
            mvFiles = [mvFilesTemp2[0]] + [fetchFolder + '/ncbi_dataset/' + tmp for tmp in mvFilesTemp2[1:]]

            for mvFile in mvFiles:
                command = f'mv "{mvFile}" "{chromosomesFolder}" -v'

                moveShell = os.popen(command)
                _ = moveShell.read()
                moveShell.close()
        except KeyboardInterrupt:
            sys.exit()
        except Exception as e:
            print(tabulation + red('ERROR:') + str(e))
            sys.exit()

        readyInfos[organismName] = fetchInfos[organismName]
        readyInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {chromosomesFolder}').read()
        readyInfos[organismName]['popular-name'] = fetchInfos[organismName]['popular-name']
        readyInfos[organismName]['kingdom'] = fetchInfos[organismName]['kingdom']

        shellTouch = os.popen('> "../rehydrated.status"')
        _ = shellTouch.read()
        shellTouch.close()

        if verbose:
            pe(green('Rehydrated'))
            print()

        os.chdir(globalGenomesPath)
    
    addToReadyFile(readyInfos, __readyFile)

    totalTemp = len(list(fetchInfos.keys()))
    rehydTemp = totalTemp - skipped - found - lost
    print(
        f'{tabulation}{numberP(rehydTemp) + magenta(f" genomes rehydrated") + " | " + numberP(found) + magenta(" found, ") + numberP(lost) +
        magenta(" lost & ") + numberP(skipped) + magenta(" skipped"):^185}'
    )
    separator()



def trnaScanSE(__readyFile=None, verbose=True, recycle=True):
    if __readyFile == None:
        global globalReadyFile
        readyFile = globalReadyFile
    else:
        readyFile = __readyFile

    with open(readyFile, 'r') as fileHandler:
        readyLines = fileHandler.readlines()
    
    readyInfos = {}
    for info in readyLines:
        splitted = info.strip('\n').split(' > ')
        popularName = splitted[3] if splitted[3] != 'None' else None

        readyInfos[splitted[1]] = {
            'accession': splitted[0],
            'tax-id': splitted[2],
            'popular-name': popularName,
            'chromosomes-folder': splitted[4],
            'kingdom': splitted[5]
        }

    found = 0
    skipped = 0
    print(f'{tabulation}{magenta(f"tRNAscan-SE analysis starting" + " | " + numberP(len(readyLines)) + magenta(" genomes")):^140}\n')
    totalIndexing = len(list(readyInfos.keys()))
    indexingSize = len(str(totalIndexing))
    for i, organismName in enumerate(readyInfos, 1):
        chromosomesFolder = readyInfos[organismName]['chromosomes-folder']
        popularName = readyInfos[organismName]['popular-name']
        accession = readyInfos[organismName]['accession']
        kingdom = readyInfos[organismName]['kingdom']
        name = f'{i:<{indexingSize}}/{totalIndexing}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')

        shellAnalysed = os.popen(f'if [ -e {globalGenomesPath}/{accession}/analysed.status ]; then echo "1"; else echo "0"; fi;')
        shellRecycled = os.popen(f'if [ -e {globalGenomesPath}/{accession}/recycled.status ]; then echo "1"; else echo "0"; fi;')

        if (verbose):
            print(f'{tabulation}{yellow(name)}:')

        if int(shellAnalysed.read()):
            found += 1
            if verbose:
                ps(red('Already analysed'))
                pe(magenta('Skipping organism'))
                print()
            continue
        
        try:
            os.chdir(chromosomesFolder)
            if verbose:
                ps(green('Chromosomes folder found'))
                pe(cyan('Attempting analysis'))
        except:
            if verbose:
                ps(red('Could not find chromosomes folder!'))
                pe(cyan('Skipping organism'))
                print()

            skipped += 1
            continue

        try:
            filesToAnalyse = list(glob.glob('*.fna'))
            filesToAnalyse.sort()

            for fileToAnalyse in filesToAnalyse:
                if verbose:
                    psT(cyan(fileToAnalyse))

                shellChromosomeAnalysed = os.popen(f'if [ -e {chromosomesFolder}/{fileToAnalyse[:-4]}.hits ]; then echo "1"; else echo "0"; fi;')
                if int(shellChromosomeAnalysed.read()):
                    if verbose:
                        pm(red('Already analysed'))

                    if recycle:
                        shellRecycled = os.popen(f'if [ -e {globalGenomesPath}/{accession}/recycled.status ]; then echo "1"; else echo "0"; fi;')
                        shellRecycled.read()
                        recycledBool = not bool(shellRecycled.read())
                        shellRecycled.close()
                        if recycledBool:
                            bigFile = 0 if fileToAnalyse == 'chromosome.fna' else 1
                            recycleFile(fileToAnalyse, big=bigFile)
                            if verbose:
                                pm(red('chromosome recycled'))

                    if verbose:
                        peT(magenta('Skipping chromosome'))

                    shellChromosomeAnalysed.close()
                    continue

                shellChromosomeAnalysed.close()

                shellTRNA = os.popen(f'tRNAscan-SE -{kingdom} {fileToAnalyse} -q --detail -D -Q -a {fileToAnalyse[:-4]}.hits')
                _ = shellTRNA.read()
                shellTRNA.close()

                if recycle:
                    recycleFile(fileToAnalyse)
                    if verbose:
                        pm(red('chromosome recycled'))

                if verbose:
                    peT(green('Finished'))
            
            if recycle:
                shellBig = os.popen(f'if [ -e {globalGenomesPath}/{accession}/fetched.status ]; then echo "1"; else echo "0"; fi;')
                bigFile = int(shellBig.read())
                recycleFile(None, recycleAll=True, big=bigFile)

                if verbose:
                    pprint(magenta('Everything recycled'))
        except KeyboardInterrupt:
            shellTRNA.close()

            shellRemoveHit = os.popen(f'rm {fileToAnalyse[:-4]}.hits')
            _ = shellRemoveHit.read()
            shellRemoveHit.close()

            sys.exit()
        except Exception as e:
            shellTRNA.close()

            shellRemoveHit = os.popen(f'rm {fileToAnalyse[:-4]}.hits')
            _ = shellRemoveHit.read()
            shellRemoveHit.close()
            
            print(tabulation + red('ERROR:') + str(e))
            sys.exit()

        shellTouch = os.popen('> "../analysed.status"')
        _ = shellTouch.read()
        shellTouch.close()

        if verbose:
            print(f'{tabulation}{green("Analysed")}')
            print()

        os.chdir(globalGenomesPath)

    print(
        f'{tabulation}{numberP(len(list(readyInfos.keys())) - found) + magenta(" genomes analysed | ") + numberP(found) +
        magenta(" already analysed & ") + numberP(skipped) + magenta(" skipped"):^180}'
    )
    separator()



def recycleFile(recycleFile, recycleAll=False, big=False):
    if recycleFile != None:
        shellRecycle = os.popen(f'rm {recycleFile}')
        _ = shellRecycle.read()
        shellRecycle.close()

    if recycleAll:
        shellTouchRecycle = os.popen('> "../recycled.status"')
        _ = shellTouchRecycle.read()
        shellTouchRecycle.close()

        shellRemoveDownloaded = os.popen('rm "../downloaded.status"')
        _ = shellRemoveDownloaded.read()
        shellRemoveDownloaded.close()

        if big:
            shellRemoveFetched = os.popen('rm "../fetched.status"')
            _ = shellRemoveFetched.read()
            shellRemoveFetched.close()



def findDetectedSeC(__readyFile=None, __detectedFile=None, verbose=True):
    if __readyFile == None:
        global globalReadyFile
        readyFile = globalReadyFile
    else:
        readyFile = __readyFile

    with open(readyFile, 'r') as fileHandler:
        readyLines = fileHandler.readlines()
    
    readyInfos = {}
    for info in readyLines:
        splitted = info.strip('\n').split(' > ')
        popularName = splitted[3] if splitted[3] != 'None' else None

        readyInfos[splitted[1]] = {
            'accession': splitted[0],
            'tax-id': splitted[2],
            'popular-name': popularName,
            'chromosomes-folder': splitted[4],
            'kingdom': splitted[5]
        }

    skipped = 0
    totalTRNA = 0
    foundPlus = 0
    foundMinus = 0
    detectedPlus = 0
    detectedInfos = {}

    print(f'{tabulation}{magenta(f"Detected tRNAs-SeC analysis starting" + " | " + numberP(len(readyLines)) + magenta(" genomes")):^140}\n')
    totalIndexing = len(list(readyInfos.keys()))
    indexingSize = len(str(totalIndexing))
    for i, organismName in enumerate(readyInfos, 1):
        chromosomesFolder = readyInfos[organismName]['chromosomes-folder']
        popularName = readyInfos[organismName]['popular-name']
        accession = readyInfos[organismName]['accession']
        kingdom = readyInfos[organismName]['kingdom']
        name = f'{i:<{indexingSize}}/{totalIndexing}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')

        shellAlreadyDetectedPlus = os.popen(f'if [ -e {globalGenomesPath}/{accession}/detected+.status ]; then echo "1"; else echo "0"; fi;')
        shellAlreadyDetectedMinus = os.popen(f'if [ -e {globalGenomesPath}/{accession}/detected-.status ]; then echo "1"; else echo "0"; fi;')
        if int(shellAlreadyDetectedPlus.read()):
            foundPlus += 1
            detectedInfos[organismName] = readyInfos[organismName]

            shellDetected = os.popen(f'cat {globalGenomesPath}/{accession}/chromosomes/*.hits | grep "SeC" -c')
            numberDetected = int(shellDetected.read())
            shellDetected.close()
            totalTRNA += numberDetected

            if verbose:
                print(f'{tabulation}{yellow(name)}:')
                ps(red('Already passed detection'))
                pe(magenta('Skipping organism'))
                print()

            continue
        elif int(shellAlreadyDetectedMinus.read()):
            foundMinus += 1

            if verbose:
                print(f'{tabulation}{yellow(name)}:')
                ps(red('Already detected'))
                pe(magenta('Skipping organism'))
                print()

            continue

        if (verbose):            
            print(f'{tabulation}{yellow(name)}:')
        
        try:
            os.chdir(chromosomesFolder)
            if verbose:
                ps(green('Chromosomes folder found'))
                pe(cyan('Attempting count'))
        except:
            if verbose:
                ps(red('Could not find chromosomes folder!'))
                pe(cyan('Skipping organism'))
                print()

            skipped += 1
            continue

        try:
            shellDetected = os.popen(f'cat *.hits | grep "SeC" -c')
            numberDetected = int(shellDetected.read())
            shellDetected.close()
        except KeyboardInterrupt:
            shellDetected.close()
            sys.exit()
        except Exception as e:
            shellDetected.close()
            print(tabulation + red('ERROR:') + str(e))
            sys.exit()


        if numberDetected > 0:
            detectedInfos[organismName] = readyInfos[organismName]
            detectedPlus += numberDetected
            totalTRNA += numberDetected

            shellTouch = os.popen('> "../detected+.status"')
            _ = shellTouch.read()
            shellTouch.close()
        else:
            shellTouch = os.popen('> "../detected-.status"')
            _ = shellTouch.read()
            shellTouch.close()

        if verbose:
            ps(green('Detection finished'))
            pe(numberP(numberDetected) + magenta(' tRNAs-SeC detected'))
            print()

        os.chdir(globalGenomesPath)

    addToDetectedFile(detectedInfos, __detectedFile)

    tempPassed = len(list(readyInfos.keys())) - foundPlus - foundMinus
    print(f'{tabulation}{numberP(tempPassed) + magenta(" genomes passed detection and ") + numberP(detectedPlus) + magenta(" tRNAs-SeC were found"):^150}')
    print(
        f'{tabulation}{numberP(foundPlus + foundMinus) + magenta(" already passed (") + numberP(foundPlus) + magenta("+ ") + numberP(foundMinus) + magenta("-) & ")
        + numberP(skipped) + magenta(" skipped | total tRNA-SeC found = ") + numberP(totalTRNA):^190}'
    )
    separator()



def taxonomyCollection(__readyFile=None, __taxonomyFile=None, verbose=True):
    global globalTaxonLevels

    if __readyFile == None:
        global globalReadyFile
        readyFile = globalReadyFile
    else:
        readyFile = __readyFile

    with open(readyFile, 'r') as fileHandler:
        readyLines = fileHandler.readlines()
    
    readyInfos = {}
    for info in readyLines:
        splitted = info.strip('\n').split(' > ')
        popularName = splitted[3] if splitted[3] != 'None' else None

        readyInfos[splitted[1]] = {
            'accession': splitted[0],
            'tax-id': splitted[2],
            'popular-name': popularName,
            'chromosomes-folder': splitted[4],
            'kingdom': splitted[5]
        }

    count = 0
    found = 0
    skipped = 0
    taxonomyInfos = {}

    print(f'{tabulation}{magenta(f"Taxonomy collection starting" + " | " + numberP(len(readyLines)) + magenta(" genomes")):^140}\n')
    totalIndexing = len(list(readyInfos.keys()))
    indexingSize = len(str(totalIndexing))
    for i, organismName in enumerate(readyInfos, 1):
        chromosomesFolder = readyInfos[organismName]['chromosomes-folder']
        popularName = readyInfos[organismName]['popular-name']
        accession = readyInfos[organismName]['accession']
        taxId = readyInfos[organismName]['tax-id']
        name = f'{i:<{indexingSize}}/{totalIndexing}. {organismName}' + cyan(f' ({popularName} - {accession}/{taxId})' if popularName != None else f' ({accession}/{taxId})')

        shellTaxonomyCollected = os.popen(f'if [ -e {globalGenomesPath}/{accession}/taxonomy.status ]; then echo "1"; else echo "0"; fi;')
        if int(shellTaxonomyCollected.read()):
            if verbose:
                print(f'{tabulation}{yellow(name)}:')
                ps(red('Already passed taxonomy collection'))
                pe(magenta('Skipping organism'))
                print()

            shellTaxonomy = os.popen(f'cat {globalGenomesPath}/{accession}/taxonomy.status')
            shellTaxonomyRead = shellTaxonomy.read().strip('\n')
            shellTaxonomy.close()

            taxonomyInfos[organismName] = readyInfos[organismName]
            taxonomyInfos[organismName]['taxonomy'] = {}
            taxonomyLines = shellTaxonomyRead.split(' > ')
            for line in taxonomyLines:
                level, taxon, taxonId = line.split('<')
                taxonomyInfos[organismName]['taxonomy'][level] = taxon
                taxonomyInfos[organismName]['taxonomy'][f'{level}-id'] = taxonId

            found += 1

            continue

        if (verbose):            
            print(f'{tabulation}{yellow(name)}:')

        try:
            os.chdir(f'{globalGenomesPath}/{accession}')
            if verbose:
                ps(green('Download folder found'))
                pe(cyan('Attempting collection'))
        except:
            if verbose:
                ps(red('Could not find download folder!'))
                pe(cyan('Skipping organism'))
                print()

            skipped += 1
            continue

        try:
            created = 0
            shellTaxonomy = os.popen(f'datasets summary taxonomy taxon "{taxId}" --as-json-lines 2>&1')

            summaryRead = shellTaxonomy.read().replace('true', "'true'")[:-1]
            shellTaxonomy.close()

            if 'does not match any existing taxids' in summaryRead:
                skipped += 1
                if verbose:
                    ps(red(f'No organism found with taxid {taxId}'))
                    pe(magenta('Skipping organism'))
                    print()
                continue

            summary = {}
            for j, summa in enumerate(summaryRead.split('\n')):
                try:
                    summary[j] = ast.literal_eval(summa)
                except Exception as e:
                    print(j, summa, e)
                    sys.exit()

            commandParts = []
            for level in globalTaxonLevels:
                if verbose:
                    ps(yellow(level))
                
                try:
                    taxonSummary = summary[0]['taxonomy']['classification'][level]
                    taxon, taxonId = taxonSummary['name'], taxonSummary['id']
                except Exception as e:
                    if level in str(e):
                        pm(red(f'Organism has no {level} (what are you trying to do?)'))
                        pe(magenta('Skipping taxonomic level!'))
                        continue
                    else:
                        shellTaxonomy.close()
                        print(tabulation + red('ERROR:') + str(e))
                        sys.exit()

                taxonomyInfos[organismName] = readyInfos[organismName]
                taxonomyInfos[organismName]['taxonomy'] = {}
                taxonomyInfos[organismName]['taxonomy'][level] = taxon
                taxonomyInfos[organismName]['taxonomy'][f'{level}-id'] = taxonId
                commandParts.append(f'{level}<{taxon}<{taxonId}')

                if verbose:
                    pe(green('Collected'))

            command = ' > '.join(commandParts)

            shellCreateTaxon = os.popen(f'echo "{command}" > taxonomy.status')
            _ = shellCreateTaxon.read()
            shellCreateTaxon.close()
            created = 1

            del taxon, taxonId

            count += 1
            if verbose:
                pprint(red(f'Taxonomy collected!'))
                print()
        except KeyboardInterrupt:
            shellTaxonomy.close()
            shellCreateTaxon.close()

            shellDeleteTaxon = os.popen(f'rm {level}.status')
            _ = shellDeleteTaxon.read()
            shellDeleteTaxon.close()

            sys.exit()
        except Exception as e:
            shellTaxonomy.close()
            shellCreateTaxon.close()

            if created:
                shellDeleteTaxon = os.popen(f'rm {level}.status')
                _ = shellDeleteTaxon.read()
                shellDeleteTaxon.close()

            print(tabulation + red('ERROR:') + str(e))
            sys.exit()
    
    addToTaxonomyFile(taxonomyInfos, __taxonomyFile)
    print(
        f'{tabulation}{magenta(f"Taxonomy collection finished | ") + numberP(count) + magenta(" genomes processed, ") + numberP(found) + magenta(" found & ")
        + numberP(skipped) + magenta(" skipped"):^175}'
    )
    separator()



def collectRSSU(__detectedFile=None, __taxonomyFile=None, __rSSUFile=None, verbose=True, rssFile='SILVA_138.2_SSURef.rnac'):
    global globalGenomesPath
    
    if __detectedFile == None:
        global globalDetectedFile
        detectedFile = globalDetectedFile
    else:
        detectedFile = __detectedFile
    
    if __taxonomyFile == None:
        global globalTaxonomyFile
        taxonomyFile = globalTaxonomyFile
    else:
        taxonomyFile = __taxonomyFile

    with open(detectedFile) as fileHandler:
        detectedLines = fileHandler.readlines()

    with open(taxonomyFile) as fileHandler:
        taxonomyLines = fileHandler.readlines()

    detectedNames = [splitted.strip('\n').split(' > ')[1] for splitted in detectedLines]

    taxonomyInfos = {}
    for line in taxonomyLines:
        accession, organismName, taxId, popularName, chrFolder, taxonomy = line.strip('\n').split(' > ')
        popularName = popularName if popularName != 'None' else None

        if not organismName in detectedNames: continue

        taxonomyInfos[organismName] = {
            'accession': accession,
            'tax-id': taxId,
            'popular-name': popularName,
            'chromosomes-folder': chrFolder
        }

        taxonomyInfos[organismName]['taxonomy'] = {}
        for info in taxonomy.split(' | '):
            level, name, taxId = info.split('<')
            taxonomyInfos[organismName]['taxonomy'][level] = {'name': name, 'tax-id': taxId}

    print(f'{tabulation}{magenta(f"rRNAs collection starting" + " | " + numberP(len(detectedNames)) + magenta(" genomes")):^140}\n')
    totalIndexing = len(list(taxonomyInfos.keys()))
    indexingSize = len(str(totalIndexing))

    try:
        global globalGenomesPath
        os.chdir(globalGenomesPath + '/../16-18S DB')
        if verbose:
            print(f'{green("16S and 18S Databases folder found!") + magenta(" | ") + magenta("Attempting collection"):^215}')
            print()
    except Exception as e:
        print(tabulation + red('ERROR:') + str(e))
        sys.exit()

    lost = 0
    found = 0
    collected = 0
    rSSUInfos = {}
    for i, organismName in enumerate(taxonomyInfos, 1):
        popularName = taxonomyInfos[organismName]['popular-name']
        accession = taxonomyInfos[organismName]['accession']
        taxId = taxonomyInfos[organismName]['tax-id']

        namePart = cyan(f' ({popularName} - {accession}/{taxId})' if popularName != None else f' ({accession}/{taxId})')
        name = f'{i:<{indexingSize}}/{totalIndexing}. {organismName}' + namePart

        if (verbose):            
            pprint(yellow(name) + ':')

        shellSSUCollectedPlus = os.popen(f'if [ -e {globalGenomesPath}/{accession}/SSUSequence+.fasta ]; then echo "1"; else echo "0"; fi;')
        shellSSUCollectedMinus = os.popen(f'if [ -e {globalGenomesPath}/{accession}/SSUSequence-.fasta ]; then echo "1"; else echo "0"; fi;')
        if int(shellSSUCollectedPlus.read()):
            if verbose:
                ps(red('SSU already collected'))
                pe(magenta('Skipping organism'))
                print()

            found += 1
            rSSUInfos[organismName] = taxonomyInfos[organismName]
            continue
        elif int(shellSSUCollectedMinus.read()):
            if verbose:
                ps(red('Already marked as not found'))
                pe(magenta('Skipping organism'))
                print()

            lost += 1
            continue

        try:
            shellCollectSSU = os.popen(f'grep "\t{taxId}\t" {rssFile}')
            readCollectRSS = shellCollectSSU.read().replace('--\n', '')
            shellCollectSSU.close()
        except KeyboardInterrupt:
            shellCollectSSU.close()
            sys.exit()
        except Exception as e:
            print(tabulation + red('ERROR:') + str(e))
            sys.exit()

        result = False
        if readCollectRSS != '':
            rawAllInfos = readCollectRSS.split('\n')[:-1]
            
            allInfos = {}
            for j, entry in enumerate(rawAllInfos):
                splitted = entry.split('\t')
                allInfos[j] = {globalRSSUFileHeader[k]: splitted[k] for k in range(len(globalRSSUFileHeader))}

            result, header, sequence = checkTaxonomyRSSU(organismName, allInfos, taxonomyInfos)
            
        if result:
            pprint(green('Found!'))
            rSSUInfos[organismName] = taxonomyInfos[organismName]

            try:
                shellAddSSU = os.popen(f'echo "{header}\n{sequence}" > {globalGenomesPath}/{accession}/SSUSequence+.fasta')
                _ = shellAddSSU.read()
                shellAddSSU.close()
            except KeyboardInterrupt:
                shellAddSSU.close()
                
                shellRemoveSSU = os.popen(f'rm {globalGenomesPath}/{accession}/SSUSequence+.fasta')
                _ = shellRemoveSSU.read()
                shellRemoveSSU.close()

                sys.exit()
            except Exception as e:
                shellAddSSU.close()
                
                shellRemoveSSU = os.popen(f'rm {globalGenomesPath}/{accession}/SSUSequence+.fasta')
                _ = shellRemoveSSU.read()
                shellRemoveSSU.close()

                print(tabulation + red('ERROR:') + str(e))
                sys.exit()

            collected += 1
        else:    
            pprint(red('Not found!'))

            try:
                shellAddSSU = os.popen(f'> {globalGenomesPath}/{accession}/SSUSequence-.fasta')
                _ = shellAddSSU.read()
                shellAddSSU.close()
            except KeyboardInterrupt:
                shellAddSSU.close()
                
                shellRemoveSSU = os.popen(f'rm {globalGenomesPath}/{accession}/SSUSequence-.fasta')
                _ = shellRemoveSSU.read()
                shellRemoveSSU.close()

                sys.exit()
            except Exception as e:
                shellAddSSU.close()
                
                shellRemoveSSU = os.popen(f'rm {globalGenomesPath}/{accession}/SSUSequence-.fasta')
                _ = shellRemoveSSU.read()
                shellRemoveSSU.close()

                print(tabulation + red('ERROR:') + str(e))
                sys.exit()

            lost += 1

        if verbose:
            print()

    addToRSSUFile(rSSUInfos, __rSSUFile)

    print(
        f'{tabulation}{magenta(f"rRNAs (rSSU) collection finished | ") + numberP(len(detectedLines)) + magenta(" genomes processed, ") + numberP(collected)
        + magenta(" collected, ") + numberP(found) + magenta(" found & ") + numberP(lost) + magenta(" lost"):^190}'
    )
    separator()



def processAndMetadata(__rSSUFile=None, __taxonomyFile=None, __processedFile=None, __processedRSSUFile=None, __metadataFile=None, verbose=True):
    global globalTRNACount, globalTaxonLevels, globalGenomesPath



    if __rSSUFile == None:
        global globalCollectedRSSUFile
        rSSUFile = globalCollectedRSSUFile
    else:
        rSSUFile = __rSSUFile

    if __taxonomyFile == None:
        global globalTaxonomyFile
        taxonomyFile = globalTaxonomyFile
    else:
        taxonomyFile = __taxonomyFile



    with open(rSSUFile, 'r') as fileHandler:
        rSSULines = fileHandler.readlines()

    with open(taxonomyFile, 'r') as fileHandler:
        taxonLines = fileHandler.readlines()
    


    rSSUInfos = {}
    for info in rSSULines:
        splitted = info.strip('\n').split(' > ')
        popularName = splitted[3] if splitted[3] != 'None' else None

        rSSUInfos[splitted[1]] = {
            'accession': splitted[0],
            'tax-id': splitted[2],
            'popular-name': popularName,
            'chromosomes-folder': splitted[4]
        }
    
    taxonomyInfos = {}
    for info in taxonLines:
        splitted = info.strip('\n').split(' > ')
        popularName = splitted[3] if splitted[3] != 'None' else None
        taxonomy = splitted[5].split(' | ')

        taxonomyInfos[splitted[1]] = {
            'accession': splitted[0],
            'tax-id': splitted[2],
            'popular-name': popularName,
            'chromosomes-folder': splitted[4]
        }
    
        taxonomyInfos[splitted[1]]['taxonomy'] = {}
        for __taxonInfo in taxonomy:
            __taxonLevel, __taxonName, __taxonId = __taxonInfo.split('<')
            taxonomyInfos[splitted[1]]['taxonomy'][__taxonLevel] = __taxonName
            taxonomyInfos[splitted[1]]['taxonomy'][f'{__taxonLevel}-id'] = __taxonId



    count = 0
    tRNAs = 0
    skipped = 0
    totalTRNAs = 0
    processedInfos = {}
    processedRSSUInfos = {}
    metadataInfos = {}



    print(f'{tabulation}{magenta(f"Processing tRNAs-SeC starting" + " | " + numberP(len(rSSULines)) + magenta(" genomes")):^140}\n')
    totalIndexing = len(list(rSSUInfos.keys()))
    indexingSize = len(str(totalIndexing))
    for i, organismName in enumerate(rSSUInfos, 1):
        chromosomesFolder = rSSUInfos[organismName]['chromosomes-folder']
        popularName = rSSUInfos[organismName]['popular-name']
        accession = rSSUInfos[organismName]['accession']
        taxId = rSSUInfos[organismName]['tax-id']
        name = f'{i:<{indexingSize}}/{totalIndexing}. {organismName}' + cyan(f' ({popularName} - {accession}/{taxId})' if popularName != None else f' ({accession}/{taxId})')

        if (verbose):            
            print(f'{tabulation}{yellow(name)}:')
        
        try:
            os.chdir(chromosomesFolder)
            if verbose:
                ps(green('Chromosomes folder found'))
                pe(cyan('Attempting extraction'))
        except:
            if verbose:
                ps(red('Could not find chromosomes folder!'))
                pe(cyan('Skipping organism'))
                print()

            skipped += 1
            continue

        try:
            shellReads = os.popen('grep "SeC" *.hits -l')
            filesToRead = shellReads.read()[:-1].split('\n')
            shellReads.close()

            hScore = 0

            for fileToRead in filesToRead:
                shellDetected = os.popen(f'grep "SeC" "{fileToRead}" -A 2')
                detectedTRNAs = shellDetected.read()[:-1].replace('--\n', '').split('>')[1:]
                shellDetected.close()

                shellRSSU = os.popen(f'cat {globalGenomesPath}/{accession}/SSUSequence+.fasta')
                rSSUSequenceInfos = shellRSSU.read().split('\n')
                shellRSSU.close()

                rName, rType = rSSUSequenceInfos[0].split(';')
                rSSUSequence = rSSUSequenceInfos[1]

                totalTRNAs += int(os.popen(f'grep "SeC" {fileToRead} -c').read().strip())

                if rSSUSequence == '':
                    if verbose:
                        ps(red('Could not find rSSU sequence!'))
                        pe(cyan('Skipping organism'))
                        print()

                    skipped += 1
                    continue

                for j, aux in enumerate(detectedTRNAs, 1):
                    part = aux.split('\n')
                    headerInfos = part[0].split(' ')
                    tRNASequence = ''.join(part[1:])

                    if tRNASequence == '':
                        continue

                    chromosomeState, chromosomeNumber, tRNANumber = headerInfos[0][1:].split('.')
                    chromosomePosition = headerInfos[1].split(':')[1]
                    strand, size, score = headerInfos[2], headerInfos[5], float(headerInfos[8])

                    if score > hScore:
                        if hScore == 0: tRNAs += 1
                        hScore = score
                        headerFinal = f'>{i}'

                        processedRSSUInfos[organismName] = {'info': rSSUInfos[organismName], 'header': headerFinal, 'sequence': rSSUSequence}
                        processedInfos[organismName] = {'info': rSSUInfos[organismName], 'header': headerFinal, 'sequence': tRNASequence}
                        metadataInfos[headerFinal[1:]] = {
                            'taxonomy': taxonomyInfos[organismName]['taxonomy'],
                            'mitochondrial': 'MT' in fileToRead,
                            'rRNA-type': rType,
                            'tRNA-number': j,
                            'tax-id': taxId,
                            'score': score
                        }
            count += 1
        except KeyboardInterrupt:
            shellDetected.close()
            sys.exit()
        except Exception as e:
            shellDetected.close()
            print(tabulation + red('ERROR:') + str(e))
            sys.exit()

        if verbose:
            print(f'{tabulation}{red("Processing finished")}')
            print()

    addToProcessedFile(processedInfos, __processedFile)
    addToProcessedRSSUFile(processedRSSUInfos, __processedRSSUFile)
    addToMetadataFile(metadataInfos, __metadataFile)

    print(
        f'{tabulation}{magenta("tRNAs-SeC processing ended | ") + numberP(count) + magenta(" genomes processed (") + numberP(tRNAs) + magenta("/") + numberP(totalTRNAs) +
        magenta(" tRNAs-SeC) & ") + numberP(skipped) + magenta(" skipped"):^205}'
    )
    separator()



def taxonAnalysisFunc(__level='all', __taxonomyFile=None, __detectedFile=None, verbose=True, debug=False, sequential=False, plot=False):
    if __taxonomyFile == None:
        global globalTaxonomyFile
        taxonomyFile = globalTaxonomyFile
    else:
        taxonomyFile = __taxonomyFile

    if __detectedFile == None:
        global globalDetectedFile
        detectedFile = globalDetectedFile
    else:
        detectedFile = __detectedFile

    with open(taxonomyFile, 'r') as fileHandler:
        taxonLines = fileHandler.readlines()

    with open(detectedFile, 'r') as fileHandler:
        detectedLines = fileHandler.readlines()
    
    taxonInfos = {}
    for info in taxonLines:
        splitted = info.strip('\n').split(' > ')
        popularName = splitted[3] if splitted[3] != 'None' else None
        taxonomy = splitted[5].split(' | ')

        taxonInfos[splitted[1]] = {
            'accession': splitted[0],
            'tax-id': splitted[2],
            'popular-name': popularName,
            'chromosomes-folder': splitted[4]
        }
    
        taxonInfos[splitted[1]]['taxonomy'] = {}
        for __taxonInfo in taxonomy:
            __taxonLevel, __taxonName, __taxonId = __taxonInfo.split('<')
            taxonInfos[splitted[1]]['taxonomy'][__taxonLevel] = __taxonName
            taxonInfos[splitted[1]]['taxonomy'][f'{__taxonLevel}-id'] = __taxonId

    detectedInfos = {}
    for info in detectedLines:
        splitted = info.strip('\n').split(' > ')
        popularName = splitted[3] if splitted[3] != 'None' else None

        detectedInfos[splitted[1]] = {
            'accession': splitted[0],
            'tax-id': splitted[2],
            'popular-name': popularName,
            'chromosomes-folder': splitted[4],
            'kingdom': splitted[5]
        }

    if __level == 'all':
        global globalTaxonLevels
        levels = globalTaxonLevels
    else:
        levels = [__l for __l in __level.split(',')]
        if ((sequential) and (len(levels) == 1)):
            global globalTaxonLevelsSequence
            levels = globalTaxonLevelsSequence[levels[0]]

    print(f'{tabulation}{magenta(f"Taxon ({__level}) analysis starting" + " | " + numberP(len(taxonLines)) + magenta(" genomes")):^140}\n')
    taxonAnalysis = defaultdict(lambda: {'total': 0, 'found': 0, 'percentage': None})
    taxonAnalysisUnique = defaultdict(lambda: {'total': 0, 'found': 0, 'percentage': None})
    taxonAnalysisUniqueCount = {'total': {}, 'found': {}}
    totalIndexing = len(list(taxonInfos.keys()))
    indexingSize = len(str(totalIndexing))
    for i, organismName in enumerate(taxonInfos, 1):
        found = organismName in detectedInfos
        foundTxt = green('Found') if found else red('Not found')
        popularName = taxonInfos[organismName]['popular-name']
        accession = taxonInfos[organismName]['accession']
        taxId = taxonInfos[organismName]['tax-id']

        namePart = cyan(f' ({popularName} - {accession}/{taxId}) {foundTxt}' if popularName != None else f' ({accession}/{taxId}) {foundTxt}')
        name = f'{i:<{indexingSize}}/{totalIndexing}. {organismName}' + namePart

        if (verbose):            
            pprint(yellow(name) + ':')
            psTaxa()

        for i, level in enumerate(levels, 1):
            if not level in taxonInfos[organismName]['taxonomy']:
                taxon = 'None'
            else:
                taxon = taxonInfos[organismName]['taxonomy'][level]
            
            if verbose:
                if i == len(levels):
                    print(cyan(taxon))
                else:
                    pmTaxa(cyan(taxon))
            
            if not taxon in taxonAnalysis[level]:
                taxonAnalysis[level][taxon] = {'total': 0, 'found': 0, 'percentage': None}
                taxonAnalysisUnique[level][taxon] = {'total': 0, 'found': 0, 'percentage': None}
                taxonAnalysisUniqueCount['total'][level] = []
                taxonAnalysisUniqueCount['found'][level] = []

            taxonAnalysis[level]['total'] += 1
            taxonAnalysis[level][taxon]['total'] += 1

            if not level in taxonAnalysisUniqueCount['total']:
                taxonAnalysisUniqueCount['total'][level] = []
            if not taxonInfos[organismName]['taxonomy']['species'] in taxonAnalysisUniqueCount['total'][level]:
                taxonAnalysisUnique[level]['total'] += 1
                taxonAnalysisUnique[level][taxon]['total'] += 1
                taxonAnalysisUniqueCount['total'][level].append(taxonInfos[organismName]['taxonomy']['species'])

            if found:
                taxonAnalysis[level]['found'] += 1
                taxonAnalysis[level][taxon]['found'] += 1

                if not level in taxonAnalysisUniqueCount['found']:
                    taxonAnalysisUniqueCount['found'][level] = []
                if not taxonInfos[organismName]['taxonomy']['species'] in taxonAnalysisUniqueCount['found'][level]:
                    taxonAnalysisUnique[level]['found'] += 1
                    taxonAnalysisUnique[level][taxon]['found'] += 1
                    taxonAnalysisUniqueCount['found'][level].append(taxonInfos[organismName]['taxonomy']['species'])

        if verbose:
            print()
            print()

    for level in levels:
        for phy in taxonAnalysis[level]:
            if phy in ['found', 'total', 'percentage', 'counted']:
                continue

            found = taxonAnalysis[level][phy]['found']
            total = taxonAnalysis[level][phy]['total']
            percentage = round(100 * (found / total), 2)
            taxonAnalysis[level][phy]['percentage'] = percentage

        found = taxonAnalysis[level]['found']
        total = taxonAnalysis[level]['total']
        percentage = round(100 * (found / total), 2)
        taxonAnalysis[level]['percentage'] = percentage

    for level in levels:
        for phy in taxonAnalysisUnique[level]:
            if phy in ['found', 'total', 'percentage', 'counted']:
                continue

            found = taxonAnalysisUnique[level][phy]['found']
            total = taxonAnalysisUnique[level][phy]['total']
            percentage = round(100 * (found / total), 2)
            taxonAnalysisUnique[level][phy]['percentage'] = percentage

        found = taxonAnalysisUnique[level]['found']
        total = taxonAnalysisUnique[level]['total']
        percentage = round(100 * (found / total), 2)
        taxonAnalysisUnique[level]['percentage'] = percentage

    if (((debug) or (__level != 'all')) and (verbose)):
        pprint(red('ORGANISMS!'))
        pretty(taxonAnalysis)

        print()
        
        pprint(red('UNIQUE (SPECIES)!'))
        pretty(taxonAnalysisUnique)

        if verbose:
            print()

    if plot:
        plotTaxonAnalysis(taxonAnalysis, levels, unique=False, show=plot)
        plotTaxonAnalysis(taxonAnalysisUnique, levels, unique=True, show=plot)

    print(
        f'{tabulation}{magenta(f"Taxon ({', '.join(levels)}) analysis finished"):^120}'
    )
    separator()



def alignMAFFT(__processedFile=None, __processedRSSUFile=None, __alignFile=None, __alignRSSUFile=None, progress=False, verbose=True):
    if __processedFile == None:
        global globalProcessedFile
        processedFile = globalProcessedFile
    else:
        processedFile = __processedFile
    
    if __processedRSSUFile == None:
        global globalProcessedRSSUFile
        processedRSSUFile = globalProcessedRSSUFile
    else:
        processedRSSUFile = __processedRSSUFile

    if __alignFile == None:
        global globalAlignFile
        alignFile = globalAlignFile
    else:
        alignFile = __alignFile

    if __alignRSSUFile == None:
        global globalAlignRSSUFile
        alignRSSUFile = globalAlignRSSUFile
    else:
        alignRSSUFile = __alignRSSUFile

    tRNACount = int(os.popen(f'grep ">" {processedFile} -c').read().strip())

    print(f'{tabulation}{magenta(f"MAFFT alignment starting" + " | " + numberP(tRNACount) + magenta(" tRNAs-SeC")):^140}\n')
    try:
        args = '--auto --reorder'
        if not progress:
            args += ' --quiet'

        if verbose:
            pprint(yellow('tRNA-SeC:'))
            ps(green('starting'))

        shellMAFFT = os.popen(f'mafft {args} "{processedFile}" > "{alignFile}"')
        _ = shellMAFFT.read()
        shellMAFFT.close()

        if verbose:
            pe(red('ended'))
            print()
            pprint(yellow('rSSU:'))
            ps(green('starting'))

        shellMAFFTRSSU = os.popen(f'mafft {args} "{processedRSSUFile}" > "{alignRSSUFile}"')
        _ = shellMAFFTRSSU.read()
        shellMAFFTRSSU.close()

        if verbose:
            pe(red('ended'))
            print()
    except KeyboardInterrupt:
        sys.exit()
    except Exception as e:
        print(tabulation + red('ERROR:') + str(e))
        sys.exit()
    
    print(f'{tabulation}{magenta(f"MAFFT alignment ended"):^120}')
    separator()