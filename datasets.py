# Imports
from collections import defaultdict, deque
from termcolor import colored
import re, os, pickle, sys
from json import dumps



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
pe = lambda x: print(f'{x}')

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

def checkSize(species, speciesName, bigfile, fileName='genome', verbose=1):
    calculated, calculatedUnit = species[speciesName]['filesize-unit']
    
    if not bigfile:
        realText = os.popen(f'du -sh "{fileName}.fasta"').read()
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
            popularName = fetchInfos[name]['popular-name']
            kingdom = fetchInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {popularName} > {fetchFolderPath} > {kingdom}\n')

def addToReadyFile(readyInfos):
    global readyFile

    if not createdReadyFile:
        createReadyFile()

    with open(readyFile, 'a') as fileHandler:
        for name in readyInfos:
            accession = readyInfos[name]['accession']
            readyFolderPath = readyInfos[name]['ready-folder'].strip('\n')
            popularName = readyInfos[name]['popular-name']
            kingdom = readyInfos[name]['kingdom']

            fileHandler.write(f'{accession} > {name} > {popularName} > {readyFolderPath} > {kingdom}\n')



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
    species = defaultdict(lambda: {'accession': None, 'level': None})

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
                species[speciesName]['level'] = report['assembly_info']['assembly_level']
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
def downloadGenomes(organisms, verbose=1, checkSizes=0, progressbar=0, sizeLimit=15, zipFile='genome', genomeFile='genome', fetchFolder='fetchFolder'):
    global genomesPath, createdGenomesPath, readyFile, createdReadyFile

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
        name = f'{i}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')
        arguments = ''

        if (checkSizes or progressbar or verbose):
            print(f'{tabulation}{yellow(name)}:')

        bigfile = (not (unit in ['B', 'KB', 'MB']) or ((size > sizeLimit) and not (unit in ['B', 'KB'])))
        shellDownloaded = os.popen(f'if [ -e {genomesPath + "/" + accession}/downloaded.status ]; then echo "1"; else echo "0"; fi;')
        if bigfile:
            shellFetched = os.popen(f'if [ -e {genomesPath + "/" + accession}/fetched.status ]; then echo "1"; else echo "0"; fi;')
            if int(shellFetched.read()):
                found += 1
                if (checkSizes or progressbar or verbose):
                    ps(green('Already fetched'))
                    pe(magenta('Skipping organism'))
                    print()

                readyInfos[organismName] = organisms[organismName]
                readyInfos[organismName]['popular-name'] = organism['popular-name']
                readyInfos[organismName]['ready-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession}').read()
                readyInfos[organismName]['kingdom'] = organism['kingdom']
                continue
            elif int(shellDownloaded.read()):
                if (checkSizes or progressbar or verbose):
                    ps(green('Already downloaded'))
                    pe(magenta('Skipping organism'))
                    print()

                fetchInfos[organismName] = organisms[organismName]
                fetchInfos[organismName]['fetch-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + fetchFolder}').read()
                fetchInfos[organismName]['popular-name'] = organism['popular-name']
                continue

            ps(red(f'The estimated filesize is bigger than {sizeLimit} MB!'))
            pe(magenta('Downloading dehydrated files'))
            arguments += ' --dehydrated'
        elif int(shellDownloaded.read()):
            found += 1
            if (checkSizes or progressbar or verbose):
                ps(green('Already downloaded'))
                pe(magenta('Skipping organism'))
                print()

            readyInfos[organismName] = organisms[organismName]
            readyInfos[organismName]['popular-name'] = organism['popular-name']
            readyInfos[organismName]['ready-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession}').read()
            readyInfos[organismName]['kingdom'] = organism['kingdom']
            continue

        if not progressbar:
            arguments += ' --no-progressbar'

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
            zipShell = os.popen(f'unzip -p "{zipFile}.zip" "ncbi_dataset/data/G*" > "{genomeFile}.fasta"')
            _ = zipShell.read()
            zipShell.close()

            readyInfos[organismName] = organisms[organismName]
            readyInfos[organismName]['popular-name'] = organism['popular-name']
            readyInfos[organismName]['ready-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession}').read()
            readyInfos[organismName]['kingdom'] = organism['kingdom']
        else:
            zipShell = os.popen(f'unzip -o "{zipFile}.zip" -d "{fetchFolder}"')
            _ = zipShell.read()
            zipShell.close()

            fetchInfos[organismName] = organisms[organismName]
            fetchInfos[organismName]['fetch-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession + "/" + fetchFolder}').read()
            fetchInfos[organismName]['popular-name'] = organism['popular-name']
        
        os.system(f'touch "downloaded.status"')
        
        if verbose:
            pe(green('Downloaded'))

        if checkSizes:
            sizes = checkSize(organisms, organismName, bigfile, fileName=genomeFile)

        if verbose:
            print()

        os.chdir(genomesPath)
        
    addToReadyFile(readyInfos)
    addToFetchFile(fetchInfos)

    totalTemp = len(list(organisms.keys()))
    fetchTemp = len(list(fetchInfos.keys()))
    downloadedTemp = totalTemp - fetchTemp - found
    print(
        f'{tabulation}{magenta(f"GENOMES DOWNLOAD FINISHED") + " | " + numberP(downloadedTemp) +
        magenta(f" GENOMES DOWNLOADED, ") + numberP(found) + magenta(" FOUND & ") + numberP(fetchTemp) + magenta(" MORE FOR REHYDRATION"):^175}'
    )
    separator()



def downloadFetch(file, verbose=1, genomeFile='genome'):
    with open(genomesPath + '/' + file + '.data', 'r') as fileHandler:
        fetchLines = fileHandler.readlines()
    
    fetchInfos = {}
    readyInfos = {}
    for info in fetchLines:
        splitted = info.strip('\n').split(' > ')
        fetchInfos[splitted[1]] = {'accession': splitted[0], 'popular-name': splitted[2], 'fetch-folder': splitted[3], 'kingdom': splitted[4]}

    found = 0
    skipped = 0
    print(f'{tabulation}{magenta(f"GENOMES REHYDRATION COMENCING"):^120}\n')
    for i, organismName in enumerate(fetchInfos, 1):
        fetchFolder = fetchInfos[organismName]['fetch-folder']
        popularName = fetchInfos[organismName]['popular-name']
        accession = fetchInfos[organismName]['accession']
        name = f'{i}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')

        shellFetched = os.popen(f'if [ -e {genomesPath + "/" + accession}/fetched.status ]; then echo "1"; else echo "0"; fi;')
        if int(shellFetched.read()):
            found += 1
            if verbose:
                print(f'{tabulation}{yellow(name)}:')
                ps(green('Already fetched'))
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
                pe(cyan('Skiping organism'))
                print()

            skipped += 1
            continue

        try:
            rehydrationShell = os.popen(f'datasets rehydrate --directory . --max-workers 30 --no-progressbar')
            _ = rehydrationShell.read()
            rehydrationShell.close()

            mvFile = fetchFolder + '/ncbi_dataset/' + os.popen(f'datasets rehydrate --directory . --max-workers 30 --no-progressbar --list').read().strip('\n')
            command = f'mv "{mvFile}" "{fetchFolder}/../{genomeFile}.fasta" -v'
            moveShell = os.popen(command)
            _ = moveShell.read()
            moveShell.close()
        except KeyboardInterrupt:
            sys.exit()
        except Exception as e:
            print(tabulation + red('ERROR:') + e)
            sys.exit()

        readyInfos[organismName] = fetchInfos[organismName]
        readyInfos[organismName]['popular-name'] = fetchInfos[organismName]['popular-name']
        readyInfos[organismName]['ready-folder'] = os.popen(f'readlink -f {genomesPath + "/" + accession}').read()
        readyInfos[organismName]['kingdom'] = fetchInfos[organismName]['kingdom']

        shellTouch = os.popen('> "../fetched.status"')
        _ = shellTouch.read()
        shellTouch.close()

        if verbose:
            pe(green('Rehydrated'))
            print()

        os.chdir(genomesPath)
    
    addToReadyFile(readyInfos)

    totalTemp = len(list(fetchInfos.keys()))
    rehydTemp = totalTemp - skipped - found
    print(f'{tabulation}{numberP(rehydTemp) + magenta(f" GENOMES REHYDRATED | ") + numberP(found) + magenta(" FOUND ") + numberP(skipped) + magenta(" SKIPPED"):^165}')
    separator()

def trnaScanSE(file, verbose=1, genomeFile='genome', hitsFile='hits'):
    with open(genomesPath + '/' + file + '.data', 'r') as fileHandler:
        readyLines = fileHandler.readlines()
    
    readyInfos = {}
    for info in readyLines:
        splitted = info.strip('\n').split(' > ')
        readyInfos[splitted[1]] = {'accession': splitted[0], 'popular-name': splitted[2], 'ready-folder': splitted[3], 'kingdom': splitted[4]}

    found = 0
    skipped = 0
    print(f'{tabulation}{magenta(f"tRNAscan-SE ANALYSIS COMENCING"):^120}\n')
    for i, organismName in enumerate(readyInfos, 1):
        readyFolder = readyInfos[organismName]['ready-folder']
        popularName = readyInfos[organismName]['popular-name']
        accession = readyInfos[organismName]['accession']
        kingdom = readyInfos[organismName]['kingdom']
        name = f'{i}. {organismName}' + cyan(f' ({popularName} - {accession})' if popularName != None else f' ({accession})')

        shellAnalysed = os.popen(f'if [ -e {genomesPath + "/" + accession}/analysed.status ]; then echo "1"; else echo "0"; fi;')
        if int(shellAnalysed.read()):
            found += 1
            if verbose:
                print(f'{tabulation}{yellow(name)}:')
                print(f'{tabulation}{green("Already analysed")} | {magenta("Skipping organism")}')
                print()
            continue

        if (verbose):            
            print(f'{tabulation}{yellow(name)}:')
        
        try:
            os.chdir(readyFolder)
            if verbose:
                ps(green("Ready folder found"))
                pm(cyan("Attempting analysis"))
        except:
            if verbose:
                ps(red("Could not find ready folder!"))
                pe(cyan("Skiping organism"))
                print()

            skipped += 1
            continue

        try:
            shellTRNA = os.popen(f'tRNAscan-SE -{kingdom} {genomeFile}.fasta -q --detail -Q -a {hitsFile}.fasta')
            _ = shellTRNA.read()
            shellTRNA.close()
        except KeyboardInterrupt:
            sys.exit()
        except Exception as e:
            print(tabulation + red('ERROR:') + e)
            sys.exit()

        shellTouch = os.popen('> "analysed.status"')
        _ = shellTouch.read()
        shellTouch.close()

        if verbose:
            pe(green("Analysed"))
            print()

        os.chdir(genomesPath)

    shellGrepNumber = os.popen(f'cat GCF_*/{hitsFile}.fasta | grep "SeC" -c')
    totalHits = int(shellGrepNumber.read())
    shellGrepNumber.close()

    print(
        f'{tabulation}{numberP(len(list(readyInfos.keys()))) + magenta(" GENOMES ANALYSED AND ") + numberP(totalHits) + magenta(" tRNA-SeCs FOUND | ") +
        numberP(found) + magenta(" ALREADY ANALYSED & ") + numberP(skipped) + magenta(" SKIPPED"):^180}'
    )
    separator()