# Imports
from collections import defaultdict, deque
from termcolor import colored
from json import dumps
import os, sys, glob



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

# Consts
globalGenomesPath = None
createdGenomesPath = False

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

def pathCheckCreation(path, returnPath=False):
    if not os.path.isdir(path):
        os.mkdir(path)
    os.chdir(path)

    if returnPath:
        return os.getcwd()
    return



# Creation of the Genomes path
def createGenomesPath(__genomesPath='Genomes/'):
    global globalGenomesPath, createdGenomesPath

    __path = os.getcwd()
    __path = pathCheckCreation(__path + '/' + __genomesPath, returnPath=True)
    createdGenomesPath = True

    globalGenomesPath = __path

def createFetchFile(__fetchFile='fetch'):
    global globalFetchFile, createdFetchFile
    
    os.system(f'> {__fetchFile}.data')
    createdFetchFile = True

    globalFetchFile = os.popen(f'readlink -f {__fetchFile}.data').read()[:-1]

def createReadyFile(__readyFile='ready'):
    global globalReadyFile, createdReadyFile
    
    os.system(f'> {__readyFile}.data')
    createdReadyFile = True

    globalReadyFile = os.popen(f'readlink -f {__readyFile}.data').read()[:-1]

def createDetectedFile(__detectedFile='detected'):
    global globalDetectedFile, createdDetectedFile
    
    os.system(f'> {__detectedFile}.data')
    createdDetectedFile = True

    globalDetectedFile = os.popen(f'readlink -f {__detectedFile}.data').read()[:-1]

def createProcessedFile(__processedFile='processed'):
    global globalProcessedFile, createdProcessedFile
    
    os.system(f'> {__processedFile}.fna')
    createdProcessedFile = True

    globalProcessedFile = os.popen(f'readlink -f {__processedFile}.fna').read()[:-1]

def createAlignFile(__alignFile='align'):
    global globalAlignFile, createdAlignFile
    
    os.system(f'> {__alignFile}.fasta')
    createdAlignFile = True

    globalAlignFile = os.popen(f'readlink -f {__alignFile}.fasta').read()[:-1]

def initiate(__genomesPath='Genomes/', __fetchFile='fetch', __readyFile='ready', __detectedFile='detected', __processedFile='processed', __alignFile='align'):
    createGenomesPath(__genomesPath)
    createFetchFile(__fetchFile)
    createReadyFile(__readyFile)
    createDetectedFile(__detectedFile)
    createProcessedFile(__processedFile)
    createAlignFile(__alignFile)

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
            fetchFolderPath = fetchInfos[name]['fetch-folder'].strip('\n')
            chromosomesFolderPath = fetchInfos[name]['chromosomes-folder'].strip('\n')
            popularName = fetchInfos[name]['popular-name']
            kingdom = fetchInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {popularName} > {fetchFolderPath} > {chromosomesFolderPath} > {kingdom}\n')

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
            chromosomesFolderPath = readyInfos[name]['chromosomes-folder'].strip('\n')
            popularName = readyInfos[name]['popular-name']
            kingdom = readyInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {popularName} > {chromosomesFolderPath} > {kingdom}\n')

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
            chromosomesFolderPath = detectedInfos[name]['chromosomes-folder'].strip('\n')
            popularName = detectedInfos[name]['popular-name']
            kingdom = detectedInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {popularName} > {chromosomesFolderPath} > {kingdom}\n')

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

            fileHandler.write(f'{header}\n{sequence}\n\n')



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
            print(pcolor(packet))




# Collection of organism data (accession and assembly level)
def collectInfo(taxons, verbose=1, debug=0):
    species = defaultdict(lambda: {'accession': None})

    separator()
    print(f'{tabulation}{magenta("Data collection starting"):^120}\n')
    for i, (taxon, (popularName, kingdom)) in enumerate(sorted(taxons.items()), 1):
        command = f'datasets summary genome taxon "{taxon}" --reference --assembly-source RefSeq 2>&1'
        shell = os.popen(command)
        summaryRead = shell.read()
        name = f'{i}. {taxon}'

        try:
            summary = eval(summaryRead.replace('true', '"true"'))
            if debug:
                pretty(summary)
        except:
            if 'but no genome data is currently available for this taxon' in summaryRead:
                printCollection(name, 0, 2, popularName=popularName)
                continue
            elif 'New version of client' in summaryRead:
                printCollection('', 0, 3)
                continue
            else:
                printCollection(name, 0, 4, popularName=popularName, summaryRead=summaryRead[:300])
                continue
        
        shell.close()

        try:
            for report in summary['reports']:
                speciesName = report['organism']['organism_name']
                species[speciesName]['accession'] = report['accession']
                species[speciesName]['filesize-unit'] = sizeOf(int(report["assembly_stats"]["total_sequence_length"]))
                species[speciesName]['popular-name'] = popularName
                species[speciesName]['kingdom'] = kingdom

            if verbose:
                numberOfOrganisms = len(summary['reports'])
                printCollection(name, numberOfOrganisms, 1, popularName=popularName)
        except:
            printCollection(name, 0, 0, popularName=popularName)

    print(f'\n{tabulation}{magenta(f"Data collection ended") + " | " + numberP(len(list(species.keys()))) + magenta(" total organisms collected"):^140}')
    separator()

    if debug:
        print('Species = ', end='')
        pretty(species)
        print()

    return species

# Downloading genomes
def downloadGenomes(
        organisms,
        __genomesPath=None,
        __readyFile=None,
        __fetchFile=None,
        referenceRange=None,
        verbose=1,
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
        organismsNamesSorted = list(organisms.keys())[referenceRange[0]:referenceRange[1]]

    print(f'{tabulation}{magenta("Genomes download starting"):^120}\n')
    for i, organismName in enumerate(organismsNamesSorted, 1):
        organism = organisms[organismName]
        accession = organism['accession']
        size, unit = organism['filesize-unit']
        popularName = organism['popular-name']
        kingdom = organism['kingdom']
        name = f'{i}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')
        arguments = ''

        shellDownloaded = os.popen(f'if [ -e {genomesPath + "/" + accession}/downloaded.status ]; then echo "1"; else echo "0"; fi;')
        shellRehydrated = os.popen(f'if [ -e {genomesPath + "/" + accession}/rehydrated.status ]; then echo "1"; else echo "0"; fi;')
        shellNoChromosome = os.popen(f'if [ -e {genomesPath + "/" + accession}/noChromosome.status ]; then echo "1"; else echo "0"; fi;')

        bigfile = (not (unit in ['B', 'KB', 'MB']) or ((size > sizeLimit) and not (unit in ['B', 'KB'])))

        if (progressbar or verbose):
            print(f'{tabulation}{yellow(name)}:')

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
                readyInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + chromosomesFolder}').read()
                readyInfos[organismName]['kingdom'] = organism['kingdom']
                continue
            elif int(shellDownloaded.read()):
                if (progressbar or verbose):
                    ps(green('Already downloaded'))
                    pe(magenta('Skipping organism'))
                    print()

                fetchInfos[organismName] = organisms[organismName]
                fetchInfos[organismName]['fetch-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + fetchFolder}').read()
                fetchInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + chromosomesFolder}').read()
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
            readyInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + chromosomesFolder}').read()
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

            shellTouch = os.popen(f'touch "{organismName.replace("/", "_")}.name"')
            _ = shellTouch.read()
            shellTouch.close()

            shellDownload = os.popen(f'datasets download genome accession {accession} --filename "{zipFile}.zip"{arguments}')
            result = shellDownload.read()
            shellDownload.close()

            if 'connection reset by peer' in result:
                pm(red('Connection reseted'))
                pe(magenta('Skipping organism'))
                skipped += 1

                continue

            if 'timeout' in result:
                pm(red('Connection timeout'))
                pe(magenta('Skipping organism'))
                skipped += 1

                continue
        except KeyboardInterrupt:
            sys.exit()
        except Exception as e:
            print(tabulation + red('ERROR:') + e)
            sys.exit()
        
        if not bigfile:
            checkPath(genomesPath + '/' + accession + '/' + chromosomesFolder)
            folder = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + chromosomesFolder}').read().strip('\n')

            zipShell = os.popen(f'unzip -p "{zipFile}.zip" "ncbi_dataset/data/G*" > "{folder}/chromosome.fna"')
            _ = zipShell.read()
            zipShell.close()

            readyInfos[organismName] = organisms[organismName]
            readyInfos[organismName]['popular-name'] = organism['popular-name']
            readyInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + chromosomesFolder}').read()
            readyInfos[organismName]['kingdom'] = organism['kingdom']
        else:
            checkPath(genomesPath + '/' + accession + '/' + chromosomesFolder)

            zipShell = os.popen(f'unzip -o "{zipFile}.zip" -d "{fetchFolder}"')
            _ = zipShell.read()
            zipShell.close()

            fetchInfos[organismName] = organisms[organismName]
            fetchInfos[organismName]['fetch-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + fetchFolder}').read()
            fetchInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + chromosomesFolder}').read()
            fetchInfos[organismName]['popular-name'] = organism['popular-name']
            fetchInfos[organismName]['kingdom'] = organism['kingdom']
        
        os.system(f'touch "downloaded.status"')
        
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
        fetchInfos[splitted[1]] = {'accession': splitted[0], 'popular-name': splitted[2], 'fetch-folder': splitted[3], 'chromosomes-folder': splitted[4], 'kingdom': splitted[5]}

    lost = 0
    found = 0
    skipped = 0
    print(f'{tabulation}{magenta(f"Genomes rehydration starting"):^120}\n')
    for i, organismName in enumerate(fetchInfos, 1):
        fetchFolder = fetchInfos[organismName]['fetch-folder']
        chromosomesFolder = fetchInfos[organismName]['chromosomes-folder']
        popularName = fetchInfos[organismName]['popular-name']
        accession = fetchInfos[organismName]['accession']
        kingdom = fetchInfos[organismName]['kingdom']
        name = f'{i}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')
        arguments = ' --directory . --max-workers 30'

        shellRehydrated = os.popen(f'if [ -e {globalGenomesPath + "/" + accession}/rehydrated.status ]; then echo "1"; else echo "0"; fi;')
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
            print(tabulation + red('ERROR:') + e)
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

def trnaScanSE(__readyFile=None, verbose=1):
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
        popularName = splitted[2] if splitted[2] != 'None' else None

        readyInfos[splitted[1]] = {'accession': splitted[0], 'popular-name': popularName, 'chromosomes-folder': splitted[3], 'kingdom': splitted[4]}

    found = 0
    skipped = 0
    print(f'{tabulation}{magenta(f"tRNAscan-SE analysis starting" + " | " + numberP(len(readyLines)) + magenta(" genomes")):^140}\n')
    for i, organismName in enumerate(readyInfos, 1):
        chromosomesFolder = readyInfos[organismName]['chromosomes-folder']
        popularName = readyInfos[organismName]['popular-name']
        accession = readyInfos[organismName]['accession']
        kingdom = readyInfos[organismName]['kingdom']
        name = f'{i}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')

        shellAnalysed = os.popen(f'if [ -e {globalGenomesPath + "/" + accession}/analysed.status ]; then echo "1"; else echo "0"; fi;')
        if int(shellAnalysed.read()):
            found += 1
            if verbose:
                print(f'{tabulation}{yellow(name)}:')
                ps(green('Already analysed'))
                pe(magenta('Skipping organism'))
                print()
            continue

        if (verbose):            
            print(f'{tabulation}{yellow(name)}:')
        
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
                psT(cyan(fileToAnalyse))

                shellChromosomeAnalysed = os.popen(f'if [ -e {chromosomesFolder}/{fileToAnalyse[:-4]}.hits ]; then echo "1"; else echo "0"; fi;')
                if int(shellChromosomeAnalysed.read()):
                    peT(green('Already analysed') + ' | ' + magenta('Skipping chromosome'))
                    continue

                shellTRNA = os.popen(f'tRNAscan-SE -{kingdom} {fileToAnalyse} -q --detail -D -Q -a {fileToAnalyse[:-4]}.hits')
                _ = shellTRNA.read()
                shellTRNA.close()

                peT(green('Finished'))
        except KeyboardInterrupt:
            shellTRNA.close()
            sys.exit()
        except Exception as e:
            shellTRNA.close()
            print(tabulation + red('ERROR:') + e)
            sys.exit()

        shellTouch = os.popen('> "../analysed.status"')
        _ = shellTouch.read()
        shellTouch.close()

        if verbose:
            print(f'{tabulation}{green("Analysed")}')
            print()

        os.chdir(globalGenomesPath)

    shellGrepNumber = os.popen(f'cat GCF_*/chromosomes/*.hits | grep "SeC" -c')
    totalHits = int(shellGrepNumber.read())
    shellGrepNumber.close()

    print(
        f'{tabulation}{numberP(len(list(readyInfos.keys())) - found) + magenta(" genomes analysed and ") + numberP(totalHits) + magenta(" tRNA-SeCs found | ") +
        numberP(found) + magenta(" already analysed & ") + numberP(skipped) + magenta(" skipped"):^180}'
    )
    separator()

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
        popularName = splitted[2] if splitted[2] != 'None' else None

        readyInfos[splitted[1]] = {'accession': splitted[0], 'popular-name': popularName, 'chromosomes-folder': splitted[3], 'kingdom': splitted[4]}

    skipped = 0
    foundPlus = 0
    foundMinus = 0
    detectedInfos = {}

    print(f'{tabulation}{magenta(f"Detected tRNAs-SeC analysis starting" + " | " + numberP(len(readyLines)) + magenta(" genomes")):^140}\n')
    for i, organismName in enumerate(readyInfos, 1):
        chromosomesFolder = readyInfos[organismName]['chromosomes-folder']
        popularName = readyInfos[organismName]['popular-name']
        accession = readyInfos[organismName]['accession']
        kingdom = readyInfos[organismName]['kingdom']
        name = f'{i}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')

        shellAlreadyDetectedPlus = os.popen(f'if [ -e {globalGenomesPath + "/" + accession}/detected+.status ]; then echo "1"; else echo "0"; fi;')
        shellAlreadyDetectedMinus = os.popen(f'if [ -e {globalGenomesPath + "/" + accession}/detected-.status ]; then echo "1"; else echo "0"; fi;')
        if int(shellAlreadyDetectedPlus.read()):
            foundPlus += 1
            detectedInfos[organismName] = readyInfos[organismName]

            if verbose:
                print(f'{tabulation}{yellow(name)}:')
                ps(green('Already passed detection'))
                pe(magenta('Skipping organism'))
                print()

            continue
        elif int(shellAlreadyDetectedMinus.read()):
            foundMinus += 1

            if verbose:
                print(f'{tabulation}{yellow(name)}:')
                ps(green('Already detected'))
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
            print(tabulation + red('ERROR:') + e)
            sys.exit()


        if numberDetected > 0:
            detectedInfos[organismName] = readyInfos[organismName]

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

    totalHits = len(list(detectedInfos.keys())) - foundPlus
    tempPassed = len(list(readyInfos.keys())) - foundPlus - foundMinus
    print(
        f'{tabulation}{numberP(tempPassed) + magenta(" genomes passed detection and ") + numberP(totalHits) + magenta(" tRNAs-SeC were found | ")
        + numberP(foundPlus + foundMinus) + magenta(" already passed (") + numberP(foundPlus) + magenta("+ ") + numberP(foundMinus) + magenta("-) & ")
        + numberP(skipped) + magenta(" skipped"):^220}'
    )
    separator()

def preprocessSeC(__detectedFile=None, __processedFile=None, verbose=1):
    if __detectedFile == None:
        global globalDetectedFile
        detectedFile = globalDetectedFile
    else:
        detectedFile = __detectedFile

    with open(detectedFile, 'r') as fileHandler:
        detectedLines = fileHandler.readlines()
    
    detectedInfos = {}
    for info in detectedLines:
        splitted = info.strip('\n').split(' > ')
        popularName = splitted[2] if splitted[2] != 'None' else None

        detectedInfos[splitted[1]] = {'accession': splitted[0], 'popular-name': popularName, 'chromosomes-folder': splitted[3], 'kingdom': splitted[4]}

    skipped = 0
    processedInfos = {}

    print(f'{tabulation}{magenta(f"Preprocessing tRNAs-SeC starting" + " | " + numberP(len(detectedLines)) + magenta(" genomes")):^140}\n')
    for i, organismName in enumerate(detectedInfos, 1):
        chromosomesFolder = detectedInfos[organismName]['chromosomes-folder']
        popularName = detectedInfos[organismName]['popular-name']
        accession = detectedInfos[organismName]['accession']
        kingdom = detectedInfos[organismName]['kingdom']
        name = f'{i}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')

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
            detectedTRNA = shellDetected.read()[:-1]
            shellDetected.close()

            parts = detectedTRNA.split('\n')
            headerInfos = parts[0].split(' ')
            tRNASequence = parts[1] + parts[2]

            chromosomeState, chromosomeNumber, tRNANumber = headerInfos[0][1:].split('.')
            chromosomePosition = headerInfos[1].split(':')[1]
            strand, size, score = headerInfos[2], headerInfos[5], headerInfos[8]

            headerFinal = f'>{organismName.replace(" ", "_")} | {kingdom}_{chromosomeState}.{chromosomeNumber}:{chromosomePosition} | {strand}_{size}_{score}'

            processedInfos[organismName] = {'info': detectedInfos[organismName], 'header': headerFinal, 'sequence': tRNASequence}
        except KeyboardInterrupt:
            shellDetected.close()
            sys.exit()
        except Exception as e:
            shellDetected.close()
            print(tabulation + red('ERROR:') + e)
            sys.exit()

        if verbose:
            print(f'{tabulation}{green("Processing finished")}')
            print()

    addToProcessedFile(processedInfos, __processedFile)

    count = len(list(processedInfos.keys()))
    print(
        f'{tabulation}{magenta("tRNAs-SeC preprocessing ended | ") + numberP(count) + magenta(" genomes processed & ") + numberP(skipped) + magenta(" skipped"):^155}'
    )
    separator()

def alignMAFFT(__processedFile=None, __alignFile=None, progress=0):
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
    
    with open(globalDetectedFile, 'r') as fileHandler:
        numberGenomes = len(fileHandler.readlines())

    print(f'{tabulation}{magenta(f"MAFFT alignment starting" + " | " + numberP(numberGenomes) + magenta(" genomes")):^140}\n')
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
        print(tabulation + red('ERROR:') + e)
        sys.exit()
    
    print(f'{tabulation}{magenta(f"MAFFT alignment ended"):^120}')
    separator()