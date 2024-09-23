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
genomesPath = None
createdGenomesPath = False

fetchFile = None
createdFetchFile = False

readyFile = None
createdReadyFile = False



def printCollection(taxon, numberOfOrganisms, status, popularName=None, summaryRead=None):
    pop = f'({popularName})' if popularName != None else ''

    if status == 0:
        print(f'{tabulation}{yellow(taxon):<40} {cyan(pop):<30} {red("Failed - No reference genome"):>40} | {numberP(numberOfOrganisms)} Organisms')
    elif status == 1:
        print(f'{tabulation}{yellow(taxon):<40} {cyan(pop):<30} {green("Collected"):>40} | {numberP(numberOfOrganisms)} Organism{"s" if numberOfOrganisms > 1 else ""}')
    elif status == 2:
        print(f'{tabulation}{yellow(taxon):<40} {cyan(pop):<30} {magenta("No sequenced genome data"):>40} |')
    elif status == 3:
        print(f'{tabulation}{yellow(taxon):<40} {cyan(pop):<30} {red("Collection error"):>40} | {red(summaryRead.replace("\n", " "))}')
    else:
        print(f'{tabulation}{red("ERROR! INVALID COLLECTION STATUS!")}')


def sizeOf(num, suffix="B"):
    for unit in ("", "K", "M", "G", "T", "P", "E", "Z"):
        if abs(num) < 1024.0:
            return [float(f'{num:3.1f}'), f'{unit}{suffix}']
        num /= 1024.0
    return [float(f'{num:3.1f}'), f'Y{suffix}']

def pathCheckCreation(path, returnPath=False):
    if not os.path.isdir(path):
        os.mkdir(path)
    os.chdir(path)

    if returnPath:
        return os.getcwd()
    return



# Creation of the Genomes path
def createGenomesPath(__genomesPath='Genomes/'):
    global genomesPath, createdGenomesPath

    __path = os.getcwd()
    __path = pathCheckCreation(__path + '/' + __genomesPath, returnPath=True)
    createdGenomesPath = True

    genomesPath = __path

def createFetchFile(__fetchFile='fetch'):
    global fetchFile, createdFetchFile
    
    os.system(f'> {__fetchFile}.data')
    createdFetchFile = True

    fetchFile = os.popen(f'readlink -f {__fetchFile}.data').read()[:-1]

def createReadyFile(__readyFile='ready'):
    global readyFile, createdReadyFile
    
    os.system(f'> {__readyFile}.data')
    createdReadyFile = True

    readyFile = os.popen(f'readlink -f {__readyFile}.data').read()[:-1]



def addToFetchFile(fetchInfos):
    global fetchFile

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

def addToReadyFile(readyInfos):
    global readyFile

    if not createdReadyFile:
        createReadyFile()

    with open(readyFile, 'a') as fileHandler:
        for name in readyInfos:
            accession = readyInfos[name]['accession']
            chromosomesFolderPath = readyInfos[name]['chromosomes-folder'].strip('\n')
            popularName = readyInfos[name]['popular-name']
            kingdom = readyInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {popularName} > {chromosomesFolderPath} > {kingdom}\n')



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
    print(f'{tabulation}{magenta("DATA COLLECTION COMENCING"):^120}\n')
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
            else:
                printCollection(name, 0, 3, popularName=popularName, summaryRead=summaryRead)
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

    print(f'\n{tabulation}{magenta(f"DATA COLLECTION ENDED") + " | " + numberP(len(list(species.keys()))) + magenta(" TOTAL ORGANISM COLLECTED"):^140}')
    separator()

    if debug:
        print('Species = ', end='')
        pretty(species)
        print()

    return species

# Downloading genomes
def downloadGenomes(organisms, verbose=1, progressbar=0, sizeLimit=15, zipFile='genome', fetchFolder='fetchFolder', chromosomesFolder='chromosomes'):
    global genomesPath, createdGenomesPath, readyFile, createdReadyFile

    lost = 0
    found = 0
    fetchInfos = {}
    readyInfos = {}

    if not createdGenomesPath:
        createGenomesPath()

    if not createdReadyFile:
        createReadyFile()

    print(f'{tabulation}{magenta("GENOMES DOWNLOAD COMENCING"):^120}\n')
    for i, organismName in enumerate(organisms, 1):
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

            shellTouch = os.popen(f'touch "{organismName}.name"')
            _ = shellTouch.read()
            shellTouch.close()

            shellDownload = os.popen(f'datasets download genome accession {accession} --filename "{zipFile}.zip"{arguments}')
            _ = shellDownload.read()
            shellDownload.close()
        except KeyboardInterrupt:
            sys.exit()
        except Exception as e:
            print(tabulation + red('ERROR:') + e)
            sys.exit()
        
        if not bigfile:
            os.mkdir(genomesPath + '/' + accession + '/' + chromosomesFolder)
            folder = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + chromosomesFolder}').read().strip('\n')

            zipShell = os.popen(f'unzip -p "{zipFile}.zip" "ncbi_dataset/data/G*" > "{folder}/chromosome.fna"')
            _ = zipShell.read()
            zipShell.close()

            readyInfos[organismName] = organisms[organismName]
            readyInfos[organismName]['popular-name'] = organism['popular-name']
            readyInfos[organismName]['chromosomes-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + chromosomesFolder}').read()
            readyInfos[organismName]['kingdom'] = organism['kingdom']
        else:
            os.mkdir(genomesPath + '/' + accession + '/' + chromosomesFolder)

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
        
    addToReadyFile(readyInfos)
    addToFetchFile(fetchInfos)

    totalTemp = len(list(organisms.keys()))
    fetchTemp = len(list(fetchInfos.keys()))
    downloadedTemp = totalTemp - fetchTemp - found - lost
    print(
        f'{tabulation}{numberP(downloadedTemp) + magenta(f" GENOMES DOWNLOADED") + " | " + numberP(found)
        + magenta(" FOUND, ") + numberP(lost) + magenta(" LOST & ") + numberP(fetchTemp) + magenta(" MORE FOR REHYDRATION"):^195}'
    )
    separator()



def downloadFetch(file, verbose=1, progressbar=0):
    with open(genomesPath + '/' + file + '.data', 'r') as fileHandler:
        fetchLines = fileHandler.readlines()
    
    fetchInfos = {}
    readyInfos = {}
    for info in fetchLines:
        splitted = info.strip('\n').split(' > ')
        fetchInfos[splitted[1]] = {'accession': splitted[0], 'popular-name': splitted[2], 'fetch-folder': splitted[3], 'chromosomes-folder': splitted[4], 'kingdom': splitted[5]}

    lost = 0
    found = 0
    skipped = 0
    print(f'{tabulation}{magenta(f"GENOMES REHYDRATION COMENCING"):^120}\n')
    for i, organismName in enumerate(fetchInfos, 1):
        fetchFolder = fetchInfos[organismName]['fetch-folder']
        chromosomesFolder = fetchInfos[organismName]['chromosomes-folder']
        popularName = fetchInfos[organismName]['popular-name']
        accession = fetchInfos[organismName]['accession']
        kingdom = fetchInfos[organismName]['kingdom']
        name = f'{i}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')
        arguments = ' --directory . --max-workers 30'

        shellRehydrated = os.popen(f'if [ -e {genomesPath + "/" + accession}/rehydrated.status ]; then echo "1"; else echo "0"; fi;')
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

        os.chdir(genomesPath)
    
    addToReadyFile(readyInfos)

    totalTemp = len(list(fetchInfos.keys()))
    rehydTemp = totalTemp - skipped - found - lost
    print(
        f'{tabulation}{numberP(rehydTemp) + magenta(f" GENOMES REHYDRATED") + " | " + numberP(found) + magenta(" FOUND, ") + numberP(lost) +
        magenta(" LOST & ") + numberP(skipped) + magenta(" SKIPPED"):^185}'
    )
    separator()

def trnaScanSE(file, verbose=1):
    with open(genomesPath + '/' + file + '.data', 'r') as fileHandler:
        readyLines = fileHandler.readlines()
    
    readyInfos = {}
    for info in readyLines:
        splitted = info.strip('\n').split(' > ')
        popularName = splitted[2] if splitted[2] != 'None' else None

        readyInfos[splitted[1]] = {'accession': splitted[0], 'popular-name': popularName, 'chromosomes-folder': splitted[3], 'kingdom': splitted[4]}

    found = 0
    skipped = 0
    print(f'{tabulation}{magenta(f"tRNAscan-SE ANALYSIS COMENCING" + " | " + numberP(len(readyLines)) + magenta(" GENOMES")):^140}\n')
    for i, organismName in enumerate(readyInfos, 1):
        chromosomesFolder = readyInfos[organismName]['chromosomes-folder']
        popularName = readyInfos[organismName]['popular-name']
        accession = readyInfos[organismName]['accession']
        kingdom = readyInfos[organismName]['kingdom']
        name = f'{i}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')

        shellAnalysed = os.popen(f'if [ -e {genomesPath + "/" + accession}/analysed.status ]; then echo "1"; else echo "0"; fi;')
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

        shellTouch = os.popen('> "analysed.status"')
        _ = shellTouch.read()
        shellTouch.close()

        if verbose:
            print(f'{tabulation}{green("Analysed")}')
            print()

        os.chdir(genomesPath)

    shellGrepNumber = os.popen(f'cat GCF_*/chromosomes/*.hits | grep "SeC" -c')
    totalHits = int(shellGrepNumber.read())
    shellGrepNumber.close()

    print(
        f'{tabulation}{numberP(len(list(readyInfos.keys()))) + magenta(" GENOMES ANALYSED AND ") + numberP(totalHits) + magenta(" tRNA-SeCs FOUND | ") +
        numberP(found) + magenta(" ALREADY ANALYSED & ") + numberP(skipped) + magenta(" SKIPPED"):^180}'
    )
    separator()