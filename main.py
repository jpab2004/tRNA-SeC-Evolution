from datasets import collectInfo, downloadGenomes, downloadFetch, trnaScanSE, findDetectedSeC, preprocessSeC, alignMAFFT, phylumDetection
from datasets import initiate, pretty
import os, argparse, sys



os.system('clear')



__taxons = {
    # # Bacteria
    # 'Pseudomonas fluorescens': (None, 'B'),
    # 'Mycoplasmoides pneumoniae': ('Pneumonia', 'B'),
    # 'Escherichia coli': ('E. coli', 'B'),
    # 'Orientia tsutsugamushi': (None, 'B'),
    # 'Clostridium botulinum': (None, 'B'),
    # 'Clostridium tetani': (None, 'B'),
    # 'Clostridium difficile': (None, 'B'),
    # 'Yersinia pestis': (None, 'B'),
    # 'Chromatium Okenii': (None, 'B'),

    # # Animais
    # 'Canis lupus familiaris': ('Cachorro', 'E'),
    # 'Rattus rattus': ('Rato preto', 'E'),
    # 'Felis catus': ('Gato', 'E'),

    # # Pássaros
    # 'Alcedo atthis': ('King Fisher', 'E'),
    # 'Cariamidae': ('Seriema', 'E'),

    # # Fungos
    # 'Trichoderma': (None, 'E'),
    # 'saccharomyces cerevisiae': ('Cerveja', 'E'),
    # 'Sarcoscyphaceae': ('Bonitinho Bah', 'E'),
    # 'Wallemia': ('Fungo Monocas', 'E'),

    # Achea
    'Archaea': ['Arqueias', 'A']
}



parser=argparse.ArgumentParser()

parser.add_argument('--genomes-path', help='Path for the genomes to be stored.')
parser.add_argument('--fetch-file', help='The name for the default fetch file information.')
parser.add_argument('--ready-file', help='The name for the default ready file information.')
parser.add_argument('--detected-file', help='The name for the default detected file information.')
parser.add_argument('--processed-file', help='The name for the default processed file information.')
parser.add_argument('--align-file', help='The name for the default alignment file.')

parser.add_argument('--taxon', help='The name for the taxon to be analysed.')
parser.add_argument('--reference-range', help='Range of organisms in the taxon to be analysed.')
parser.add_argument('--range-step', help='The step of the range.')

parser.add_argument('--re-download', help='Force re download of genomes.', action='store_true')
parser.add_argument('--refseq', help='Use RefSeq expanded search.', action='store_true')

parser.add_argument('--quiet', help='Quiet printing.', action='store_false')
parser.add_argument('-QD', '--q-download', help='Quiet for download.', action='store_false')
parser.add_argument('-QF', '--q-fetch', help='Quiet for fetch.', action='store_false')
parser.add_argument('-QS', '--q-scan', help='Quiet for scan.', action='store_false')
parser.add_argument('-QDE', '--q-detected', help='Quiet for detected tRNAs-SeC.', action='store_false')
parser.add_argument('-QP', '--q-preprocess', help='Quiet for preprocess tRNAs-SeC.', action='store_false')
parser.add_argument('-QPh', '--q-phylum', help='Quiet for phylum detection.', action='store_false')

args=parser.parse_args()

genomesPath = args.genomes_path if args.genomes_path != None else 'GenomesArchaeaExpanded/'
fetchFile = args.fetch_file if args.fetch_file != None else 'fetch'
readyFile = args.ready_file if args.ready_file != None else 'ready'
detectedFile = args.detected_file if args.detected_file != None else 'detected'
processedFile = args.processed_file if args.processed_file != None else 'processed'
alignFile = args.align_file if args.align_file != None else 'align'

taxons = eval(args.taxon) if args.taxon != None else __taxons
referenceRange = int(args.reference_range) if args.reference_range != None else None
rangeStep = int(args.range_step) if args.range_step != None else 1

reDownload = args.re_download
refseq = args.refseq

quiet = args.quiet
if not quiet:
    verboseDownload = 0
    verboseFetch = 0
    verboseScan = 0
    verboseDetected = 0
    verbosePreprocess = 0
    verbosePhylum = 0
else:
    verboseDownload = args.q_download
    verboseFetch = args.q_fetch
    verboseScan = args.q_scan
    verboseDetected = args.q_detected
    verbosePreprocess = args.q_preprocess
    verbosePhylum = args.q_phylum



suffix = ''
if referenceRange != None:
    suffix = f'_{referenceRange}*{rangeStep}'

files = {
    '__genomesPath': genomesPath,
    '__fetchFile': fetchFile + suffix,
    '__readyFile': readyFile + suffix,
    '__detectedFile': detectedFile + suffix,
    '__processedFile': processedFile + suffix,
    '__alignFile': alignFile
}
initiate(**files)

species = collectInfo(taxons, archaea=refseq)
downloadGenomes(species, sizeLimit=20, referenceRange=referenceRange, rangeStep=rangeStep, reDownload=reDownload, verbose=verboseDownload)
downloadFetch(verbose=verboseFetch)
trnaScanSE(verbose=verboseScan)
findDetectedSeC(verbose=verboseDetected)
preprocessSeC(verbose=verbosePreprocess)
phylumDetection(verbose=verbosePhylum)

if referenceRange == None:
    alignMAFFT()