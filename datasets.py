# Imports
from collections import defaultdict, deque
from termcolor import colored
import os, sys, glob, ast
from json import dumps
import random, pickle



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



# Consts
globalGenomesPath = None
createdGenomesPath = False

globalSpeciesFile = 'species.pickle'

globalFetchFile = None
createdFetchFile = False

globalReadyFile = None
createdReadyFile = False

globalDetectedFile = None
createdDetectedFile = False

globalProcessedFile = None
createdProcessedFile = False

globalAlignFile = None
createdAlignFile = False

globalTaxonomyFile = None
createdTaxonomyFile = False

globalMetadataFile = None
createdMetadataFile = False

globalTRNACount = None
globalMetadataColors = None

globalTaxonLevels = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
globalTaxonLevelsCheat = {
    'superkingdom': 'k',
    'kingdom': None,
    'phylum': 'p',
    'class': 'c',
    'order': 'o',
    'family': 'f',
    'genus': 'g',
    'species': 's'
}
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



def randomColor(num=1):
    alpha = '0123456789ABCDEF'
    if num > 1:
        return [('#' + ''.join([random.choice(alpha) for _ in range(6)])) for _ in range(num)]
    return '#' + ''.join([random.choice(alpha) for _ in range(6)])



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

def createProcessedFile(__processedFile='processed', suppress=False):
    global globalProcessedFile, createdProcessedFile
    
    if not suppress: os.system(f'> {__processedFile}.fna')
    createdProcessedFile = True

    globalProcessedFile = os.popen(f'readlink -f {__processedFile}.fna').read()[:-1]

def createAlignFile(__alignFile='align', suppress=False):
    global globalAlignFile, createdAlignFile
    
    if not suppress: os.system(f'> {__alignFile}.fasta')
    createdAlignFile = True

    globalAlignFile = os.popen(f'readlink -f {__alignFile}.fasta').read()[:-1]

def createTaxonomyFile(__taxonomyFile='taxonomy', suppress=False):
    global globalTaxonomyFile, createdTaxonomyFile
    
    if not suppress: os.system(f'> {__taxonomyFile}.data')
    createdTaxonomyFile = True

    globalTaxonomyFile = os.popen(f'readlink -f {__taxonomyFile}.data').read()[:-1]

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
        __processedFile='processed',
        __alignFile='align',
        __taxonomyFile='taxonomy',
        __metadataFile='metadata',
        
        suppressDownload=False,
        suppressFetch=False,
        suppressDetected=False,
        suppressPreprocess=False,
        suppressTaxonomy=False,
        suppressAlign=False
    ):
    createGenomesPath(__genomesPath, suppress=suppressDownload)
    createFetchFile(__fetchFile, suppress=suppressFetch)
    createReadyFile(__readyFile, suppress=(suppressDownload and suppressFetch))
    createDetectedFile(__detectedFile, suppress=suppressDetected)
    createProcessedFile(__processedFile, suppress=suppressPreprocess)
    createAlignFile(__alignFile, suppress=suppressAlign)
    createTaxonomyFile(__taxonomyFile, suppress=suppressTaxonomy)
    createMetadataFile(__metadataFile, suppress=suppressPreprocess)

    return



def addToFetchFile(fetchInfos, __fetchFile=None):
    if __fetchFile == None:
        global globalFetchFile
        fetchFile = globalFetchFile
    else:
        fetchFile = __fetchFile

    if not createdFetchFile:
        createFetchFile()

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
    if __readyFile == None:
        global globalReadyFile
        readyFile = globalReadyFile
    else:
        readyFile = __readyFile

    if not createdReadyFile:
        createReadyFile()

    with open(readyFile, 'a') as fileHandler:
        for name in readyInfos:
            accession = readyInfos[name]['accession']
            taxId = readyInfos[name]['tax-id']
            chromosomesFolderPath = readyInfos[name]['chromosomes-folder'].strip('\n')
            popularName = readyInfos[name]['popular-name']
            kingdom = readyInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {taxId} > {popularName} > {chromosomesFolderPath} > {kingdom}\n')

def addToDetectedFile(detectedInfos, __detectedFile=None):
    if __detectedFile == None:
        global globalDetectedFile
        detectedFile = globalDetectedFile
    else:
        detectedFile = __detectedFile

    if not createdDetectedFile:
        createDetectedFile()

    with open(detectedFile, 'a') as fileHandler:
        for name in detectedInfos:
            accession = detectedInfos[name]['accession']
            taxId = detectedInfos[name]['tax-id']
            chromosomesFolderPath = detectedInfos[name]['chromosomes-folder'].strip('\n')
            popularName = detectedInfos[name]['popular-name']
            kingdom = detectedInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {taxId} > {popularName} > {chromosomesFolderPath} > {kingdom}\n')

def addToProcessedFile(processedInfos, __processedFile=None):
    if __processedFile == None:
        global globalProcessedFile
        processedFile = globalProcessedFile
    else:
        processedFile = __processedFile

    if not createdProcessedFile:
        createProcessedFile()

    with open(processedFile, 'a') as fileHandler:
        for name in processedInfos:
            header = processedInfos[name]['header']
            sequence = processedInfos[name]['sequence']

            fileHandler.write(f'{header}\n{sequence}\n')

def addToTaxonomyFile(taxonInfos, __taxonomyFile=None):
    global globalTaxonLevels

    if __taxonomyFile == None:
        global globalTaxonomyFile
        taxonomyFile = globalTaxonomyFile
    else:
        taxonomyFile = __taxonomyFile

    if not createdTaxonomyFile:
        createTaxonomyFile()

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
    if __metadataFile == None:
        global globalMetadataFile
        metadataFile = globalMetadataFile
    else:
        metadataFile = __metadataFile

    if not createdMetadataFile:
        createMetadataFile()

    with open(metadataFile, 'a') as fileHandler:
        fileHandler.write('ID;taxID;tRNA-SeC number\n')
        for identification in metadataInfos:
            taxId = metadataInfos[identification]['tax-id']
            taxonomy = metadataInfos[identification]['taxonomy']
            tRNANumber = metadataInfos[identification]['trna-number']

            fileHandler.write(f'{identification};{taxId};{tRNANumber}\n')



# Printing dictionaries
def pretty(value, sort_keys=True, indent=4, colorOrder=[green, blue, red, cyan, yellow, magenta]):
    color = deque(colorOrder)
    
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
def collectInfo(taxons, verbose=1, archaea=0, save=0, read=0, __speciesFile=None):
    species = defaultdict(lambda: {'accession': None})

    if __speciesFile == None:
        global globalSpeciesFile
        speciesFile = globalSpeciesFile
    else:
        speciesFile = __speciesFile

    separator()
    print(f'{tabulation}{magenta("Data collection starting"):^120}\n')

    if read:
        print(f'{tabulation}{green("Read flag detected") + " | " + magenta("Trying to read file"):^130}\n')

        with open(speciesFile, 'rb') as fileHandler:
            species = pickle.load(fileHandler)

        print(f'{tabulation}{green("Read successful"):^120}')
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
        verbose=1,
        reDownload=0,
        progressbar=0,
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



def downloadFetch(__fetchFile=None, __readyFile=None, verbose=1, progressbar=0):
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



def trnaScanSE(__readyFile=None, verbose=1, recycle=1):
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
                recycleFile(None, recycleAll=1, big=bigFile)

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

    shellGrepNumber = os.popen(f'cat {globalGenomesPath}/GCF_*/chromosomes/*.hits | grep "SeC" -c')
    totalHits = int(shellGrepNumber.read())
    shellGrepNumber.close()

    print(
        f'{tabulation}{numberP(len(list(readyInfos.keys())) - found) + magenta(" genomes analysed and ") + numberP(totalHits) + magenta(" tRNA-SeCs found | ") +
        numberP(found) + magenta(" already analysed & ") + numberP(skipped) + magenta(" skipped"):^180}'
    )
    separator()



def recycleFile(recycleFile, recycleAll=0, big=0):
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



def findDetectedSeC(__readyFile=None, __detectedFile=None, verbose=1):
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
    print(
        f'{tabulation}{numberP(tempPassed) + magenta(" genomes passed detection and ") + numberP(detectedPlus) + magenta(" tRNAs-SeC were found | ")
        + numberP(foundPlus + foundMinus) + magenta(" already passed (") + numberP(foundPlus) + magenta("+ ") + numberP(foundMinus) + magenta("-) & ")
        + numberP(skipped) + magenta(" skipped"):^220}'
    )
    separator()



def taxonCollection(__readyFile=None, __taxonomyFile=None, verbose=1):
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
    taxonInfos = {}

    print(f'{tabulation}{magenta(f"Taxonomy collection starting" + " | " + numberP(len(readyLines)) + magenta(" genomes")):^140}\n')
    totalIndexing = len(list(readyInfos.keys()))
    indexingSize = len(str(totalIndexing))
    for i, organismName in enumerate(readyInfos, 1):
        chromosomesFolder = readyInfos[organismName]['chromosomes-folder']
        popularName = readyInfos[organismName]['popular-name']
        accession = readyInfos[organismName]['accession']
        taxId = readyInfos[organismName]['tax-id']
        name = f'{i:<{indexingSize}}/{totalIndexing}. {organismName}' + cyan(f' ({popularName} - {accession}/{taxId})' if popularName != None else f' ({accession}/{taxId})')

        shellTaxonCollected = os.popen(f'if [ -e {globalGenomesPath}/{accession}/taxonomy.status ]; then echo "1"; else echo "0"; fi;')
        if int(shellTaxonCollected.read()):
            if verbose:
                print(f'{tabulation}{yellow(name)}:')
                ps(red('Already passed taxonomy collection'))
                pe(magenta('Skipping organism'))
                print()

            shellTaxon = os.popen(f'cat {globalGenomesPath}/{accession}/taxonomy.status')
            shellTaxonRead = shellTaxon.read().strip('\n')
            shellTaxon.close()

            taxonInfos[organismName] = readyInfos[organismName]
            taxonInfos[organismName]['taxonomy'] = {}
            taxonomyLines = shellTaxonRead.split(' > ')
            for line in taxonomyLines:
                level, taxon, taxonId = line.split('<')
                taxonInfos[organismName]['taxonomy'][level] = taxon
                taxonInfos[organismName]['taxonomy'][f'{level}-id'] = taxonId

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
            shellTaxon = os.popen(f'datasets summary taxonomy taxon "{taxId}" --as-json-lines 2>&1')

            summaryRead = shellTaxon.read().replace('true', "'true'")[:-1]
            shellTaxon.close()

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
                        shellTaxon.close()
                        print(tabulation + red('ERROR:') + str(e))
                        sys.exit()

                taxonInfos[organismName] = readyInfos[organismName]
                taxonInfos[organismName]['taxonomy'] = {}
                taxonInfos[organismName]['taxonomy'][level] = taxon
                taxonInfos[organismName]['taxonomy'][f'{level}-id'] = taxonId
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
            shellTaxon.close()
            shellCreateTaxon.close()

            shellDeleteTaxon = os.popen(f'rm {level}.status')
            _ = shellDeleteTaxon.read()
            shellDeleteTaxon.close()

            sys.exit()
        except Exception as e:
            shellTaxon.close()
            shellCreateTaxon.close()

            if created:
                shellDeleteTaxon = os.popen(f'rm {level}.status')
                _ = shellDeleteTaxon.read()
                shellDeleteTaxon.close()

            print(tabulation + red('ERROR:') + str(e))
            sys.exit()
    
    addToTaxonomyFile(taxonInfos, __taxonomyFile)
    print(
        f'{tabulation}{magenta(f"Taxonomy collection finished | ") + numberP(count) + magenta(" genomes processed, ") + numberP(found) + magenta(" found & ")
        + numberP(skipped) + magenta(" skipped"):^175}'
    )
    separator()



def preprocessAndMetadata(__detectedFile=None, __taxonomyFile=None, __processedFile=None, __metadataFile=None, verbose=1, debug=1):
    global globalTRNACount, globalTaxonLevels



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



    with open(detectedFile, 'r') as fileHandler:
        detectedLines = fileHandler.readlines()

    with open(taxonomyFile, 'r') as fileHandler:
        taxonLines = fileHandler.readlines()
    


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
    processedInfos = {}
    metadataInfos = {}



    print(f'{tabulation}{magenta(f"Preprocessing tRNAs-SeC starting" + " | " + numberP(len(detectedLines)) + magenta(" genomes")):^140}\n')
    totalIndexing = len(list(detectedInfos.keys()))
    indexingSize = len(str(totalIndexing))
    for i, organismName in enumerate(detectedInfos, 1):
        chromosomesFolder = detectedInfos[organismName]['chromosomes-folder']
        popularName = detectedInfos[organismName]['popular-name']
        accession = detectedInfos[organismName]['accession']
        taxId = detectedInfos[organismName]['tax-id']
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
            shellDetected = os.popen(f'cat *.hits | grep "SeC" -A 2')
            detectedTRNAs = shellDetected.read()[:-1].replace('--\n', '').split('>')[1:]
            shellDetected.close()

            if debug:
                print(detectedTRNAs, end='\n\n')

            for j, aux in enumerate(detectedTRNAs, 1):
                part = aux.split('\n')
                headerInfos = part[0].split(' ')
                tRNASequence = ''.join(part[1:])

                if tRNASequence == '':
                    continue

                if debug:
                    print(headerInfos, tRNASequence, sep='\n')

                chromosomeState, chromosomeNumber, tRNANumber = headerInfos[0][1:].split('.')
                chromosomePosition = headerInfos[1].split(':')[1]
                strand, size, score = headerInfos[2], headerInfos[5], headerInfos[8]

                # headerFinal = f'''
                #     >{taxon}.{organismName.replace(" ", "_")}.{tRNANumber} | 
                #     {taxon}_{chromosomeState}.{chromosomeNumber}:{chromosomePosition} | 
                #     {strand}_{size}_{score}
                # '''.replace('\n                    ', '').replace('\n', '')
                # headerFinal = f'>{taxon}.SeC-{j}.{organismName.replace(".", " ").replace(" ", "_")}'
                # headerFinal = f'>{taxId}.SeC-{j}.{organismName.replace(" ", "_")}'
                # headerFinal = f'>{taxon}.{organismName.replace(" ", "_")}.{taxId}.{j}'
                # headerFinal = f'>{taxId}.{j}'
                headerFinal = f'>{i}'

                metadataInfos[headerFinal[1:]] = {'taxonomy': taxonomyInfos[organismName]['taxonomy'], 'tax-id': taxId, 'trna-number': j}
                processedInfos[f'{organismName}.{j}'] = {'info': detectedInfos[organismName], 'header': headerFinal, 'sequence': tRNASequence}

                tRNAs += 1
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
    addToMetadataFile(metadataInfos, __metadataFile)
    globalTRNACount = tRNAs

    print(
        f'{tabulation}{magenta("tRNAs-SeC preprocessing ended | ") + numberP(count) + magenta(" genomes processed (") + numberP(tRNAs) + magenta(" tRNAs-SeC) & ") +
        numberP(skipped) + magenta(" skipped"):^175}'
    )
    separator()



def taxonAnalysis(__level='all', __taxonomyFile=None, __detectedFile=None, verbose=1, debug=0, sequential=False):
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
            if not level in taxonInfos[organismName]['taxonomy']: continue
            taxon = taxonInfos[organismName]['taxonomy'][level]
            if verbose:
                if i == len(levels):
                    print(cyan(taxon))
                else:
                    pmTaxa(cyan(taxon))
            
            if not taxon in taxonAnalysis[level]:
                taxonAnalysis[level][taxon] = {'total': 0, 'found': 0, 'percentage': None}

            taxonAnalysis[level]['total'] += 1
            taxonAnalysis[level][taxon]['total'] += 1
            if found:
                taxonAnalysis[level]['found'] += 1
                taxonAnalysis[level][taxon]['found'] += 1

        if verbose:
            print()
            print()

    for level in levels:
        for phy in taxonAnalysis[level]:
            if phy in ['found', 'total', 'percentage']:
                continue

            found = taxonAnalysis[level][phy]['found']
            total = taxonAnalysis[level][phy]['total']
            percentage = f'{100 * (found / total):.2f}%'
            taxonAnalysis[level][phy]['percentage'] = percentage

        found = taxonAnalysis[level]['found']
        total = taxonAnalysis[level]['total']
        percentage = f'{100 * (found / total):.2f}%'
        taxonAnalysis[level]['percentage'] = percentage

    if debug:
        pretty(taxonAnalysis)
        if verbose:
            print()

    print(
        f'{tabulation}{magenta(f"Taxon ({__level}) analysis finished"):^120}'
    )
    separator()



def collectRRS(__detectedFile=None, __taxonomyFile=None, verbose=1, king=None):
    global globalTaxonLevelsCheat

    if ((king == 'Archaea') or (king == 'Eubacteria')):
        # letterFunc = str.upper
        rssFile = 'MIMt-16S_M2c_24_4.fna'
    elif (king == 'Eukaryota'):
        # letterFunc = str.lower
        rssFile = 'General_EUK_SSU_v1.9.fasta'
    elif (king == None):
        rssFile = 'SILVA_138.2_SSURef.rnac'
    else:
        pprint(red('ERROR! INVALID KINGDOM!'))
        pprint(yellow(king))
        sys.exit()

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

    print(f'{tabulation}{magenta(f"rRNAs collection starting" + " | " + numberP(len(detectedNames)) + magenta(" genomes")):^140}')
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

    found = 0
    tot = 0
    for i, organismName in enumerate(taxonomyInfos, 1):
        popularName = taxonomyInfos[organismName]['popular-name']
        accession = taxonomyInfos[organismName]['accession']
        taxId = taxonomyInfos[organismName]['tax-id']

        namePart = cyan(f' ({popularName} - {accession}/{taxId})' if popularName != None else f' ({accession}/{taxId})')
        name = f'{i:<{indexingSize}}/{totalIndexing}. {organismName}' + namePart

        if (verbose):            
            pprint(yellow(name) + ':')

        # commandPart = ''
        # for i, level in enumerate(taxonomyInfos[organismName]['taxonomy'], 1):
        #     if level == 'kingdom': continue
        #     levelName = taxonomyInfos[organismName]['taxonomy'][level]['name']
        #     levelTaxId = taxonomyInfos[organismName]['taxonomy'][level]['tax-id']
        #     commandPart += letterFunc(globalTaxonLevelsCheat[level]) + (f'__{levelName}; ' if i<len(taxonomyInfos[organismName]['taxonomy']) else f'__{levelName}')
        # command = f'grep "{commandPart}" {rssFile}'
        command = f'grep "{taxonomyInfos[organismName]["taxonomy"]["species"]["name"]}" {rssFile} -A 1'
        # command = f'grep "{organismName}" {rssFile} -A 1'

        try:
            shellCollectRSS = os.popen(command)
            readCollectRSS = shellCollectRSS.read().replace('--\n', '')
            shellCollectRSS.close()
        except Exception as e:
            print(tabulation + red('ERROR:') + str(e))
            sys.exit()

        print(tabulation + command)
        print(tabulation + readCollectRSS)
        found += 1 if readCollectRSS != '' else 0

        if verbose:
            print()

    print(tabulation + red(found))

    print(
        f'{tabulation}{magenta(f"rRNAs collection finished"):^120}'
    )
    separator()



def alignMAFFT(__processedFile=None, __alignFile=None, progress=0):
    global globalTRNACount

    if globalTRNACount == None: tRNACount = red('UNDEFINED!')
    else: tRNACount = numberP(globalTRNACount)

    if __processedFile == None:
        global globalProcessedFile
        processedFile = globalProcessedFile
    else:
        processedFile = __processedFile

    if __alignFile == None:
        global globalAlignFile
        alignFile = globalAlignFile
    else:
        alignFile = __alignFile

    print(f'{tabulation}{magenta(f"MAFFT alignment starting" + " | " + tRNACount + magenta(" tRNAs-SeC")):^140}\n')
    try:
        args = '--auto --reorder'
        if not progress:
            args += ' --quiet'

        shellMAFFT = os.popen(f'mafft {args} "{processedFile}" > "{alignFile}"')
        _ = shellMAFFT.read()
        shellMAFFT.close()
    except KeyboardInterrupt:
        sys.exit()
    except Exception as e:
        print(tabulation + red('ERROR:') + str(e))
        sys.exit()
    
    print(f'{tabulation}{magenta(f"MAFFT alignment ended"):^120}')
    separator()