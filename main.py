from datasets import collectInfo, downloadGenomes, downloadFetch, trnaScanSE, findDetectedSeC, preprocessSeC, alignMAFFT, taxonDetection, taxonAnalysis
from datasets import initiate, pretty
import os, argparse, sys



os.system('clear')



__taxonNames = {
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
parser.add_argument('--taxon-file', help='The name for the default taxon file.')

parser.add_argument('--taxonNames', help='The names of the taxons to be analysed.')
parser.add_argument('--reference-range', help='Range of organisms in the taxon to be analysed.')
parser.add_argument('--range-step', help='The step of the range.')

parser.add_argument('--re-download', help='Force re download of genomes.', action='store_true')
parser.add_argument('--refseq', help='Use RefSeq expanded search.', action='store_true')

parser.add_argument('--quiet', help='Quiet printing.', action='store_false')
parser.add_argument('-QD', '--q-download', help='Quiet for download.', action='store_false')
parser.add_argument('-QF', '--q-fetch', help='Quiet for fetch.', action='store_false')
parser.add_argument('-QS', '--q-scan', help='Quiet for scan.', action='store_false')
parser.add_argument('-QDe', '--q-detected', help='Quiet for detected tRNAs-SeC.', action='store_false')
parser.add_argument('-QP', '--q-preprocess', help='Quiet for preprocess tRNAs-SeC.', action='store_false')
parser.add_argument('-QPh', '--q-taxon', help='Quiet for taxon detection.', action='store_false')
parser.add_argument('-QPA', '--q-taxana', help='Quiet for taxon analysis.', action='store_false')

parser.add_argument('-SD', '--suppress-download', help='Suppress download process.', action='store_false')
parser.add_argument('-SF', '--suppress-fetch', help='Suppress fetch process.', action='store_false')
parser.add_argument('-SS', '--suppress-scan', help='Suppress scan process.', action='store_false')
parser.add_argument('-SDe', '--suppress-detected', help='Suppress detected process.', action='store_false')
parser.add_argument('-SP', '--suppress-preprocess', help='Suppress preprocess.', action='store_false')
parser.add_argument('-SPh', '--suppress-taxon', help='Suppress taxon process.', action='store_false')
parser.add_argument('-SPA', '--suppress-taxana', help='Suppress taxon analysis process.', action='store_false')

args=parser.parse_args()



genomesPath = args.genomes_path if args.genomes_path != None else 'GenomesArchaeaExpanded/'
fetchFile = args.fetch_file if args.fetch_file != None else 'fetch'
readyFile = args.ready_file if args.ready_file != None else 'ready'
detectedFile = args.detected_file if args.detected_file != None else 'detected'
processedFile = args.processed_file if args.processed_file != None else 'processed'
alignFile = args.align_file if args.align_file != None else 'align'
taxonFile = args.taxon_file if args.taxon_file != None else 'phylum'

taxonNames = eval(args.taxonNames) if args.taxonNames != None else __taxonNames
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
    verboseTaxon = 0
    verboseTaxonAnalysis = 0
else:
    verboseDownload = args.q_download
    verboseFetch = args.q_fetch
    verboseScan = args.q_scan
    verboseDetected = args.q_detected
    verbosePreprocess = args.q_preprocess
    verboseTaxon = args.q_taxon
    verboseTaxonAnalysis = args.q_taxana

suppressDownload = args.suppress_download
suppressFetch = args.suppress_fetch
suppressScan = args.suppress_scan
suppressDetected = args.suppress_detected
suppressPreprocess = args.suppress_preprocess
suppressTaxon = args.suppress_taxon
suppressTaxonAnalysis = args.suppress_taxana



suffix = ''
if referenceRange != None:
    suffix = f'_{referenceRange}*{rangeStep}'

files = {
    '__genomesPath': genomesPath,
    '__fetchFile': fetchFile + suffix,
    '__readyFile': readyFile + suffix,
    '__detectedFile': detectedFile + suffix,
    '__processedFile': processedFile + suffix,
    '__alignFile': alignFile,
    '__taxonFile': taxonFile,

    'suppressDownload': not suppressDownload,
    'suppressFetch': not suppressFetch,
    'suppressDetected': not suppressDetected,
    'suppressPreprocess': not suppressPreprocess,
    'suppressTaxon': not suppressTaxon
}
initiate(**files)

species = collectInfo(taxonNames, archaea=refseq)

if suppressDownload: downloadGenomes(species, sizeLimit=20, referenceRange=referenceRange, rangeStep=rangeStep, reDownload=reDownload, verbose=verboseDownload)
if suppressFetch: downloadFetch(verbose=verboseFetch)
if suppressScan: trnaScanSE(verbose=verboseScan)
if suppressDetected: findDetectedSeC(verbose=verboseDetected)
if suppressPreprocess: preprocessSeC(verbose=verbosePreprocess)
if suppressTaxon: taxonDetection(verbose=verboseTaxon)
if suppressTaxonAnalysis: taxonAnalysis(verbose=verboseTaxonAnalysis, debug=1)

if referenceRange == None:
    alignMAFFT()