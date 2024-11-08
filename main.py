from datasets import collectInfo, downloadGenomes, downloadFetch, trnaScanSE, findDetectedSeC, preprocessAndMetadata, alignMAFFT, taxonCollection, taxonAnalysis, collectRRS
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
    # 'Archaea': ['Arqueias', 'A']
}



parser=argparse.ArgumentParser()

parser.add_argument('--genomes-path', help='Path for the genomes to be stored.')
parser.add_argument('--species-file', help='The name for the default species file information.')
parser.add_argument('--fetch-file', help='The name for the default fetch file information.')
parser.add_argument('--ready-file', help='The name for the default ready file information.')
parser.add_argument('--detected-file', help='The name for the default detected file information.')
parser.add_argument('--processed-file', help='The name for the default processed file information.')
parser.add_argument('--align-file', help='The name for the default alignment file.')
parser.add_argument('--taxonomy-file', help='The name for the default taxonomy file.')
parser.add_argument('--metadata-file', help='The name for the default metadata file.')

parser.add_argument('--read-species', help='Read species file.', action='store_true')
parser.add_argument('--save-species', help='Save species file.', action='store_true')

parser.add_argument('--taxon-names', help='The names of the taxons to be analysed.')
parser.add_argument('--taxon-level', help='The level of the taxon to be analysed.')
parser.add_argument('--sequential', help='Define if taxon analysis should consider the whole taxonomy sequence.', action='store_true')

parser.add_argument('--reference-range', help='Range of organisms in the taxon to be analysed.')
parser.add_argument('--range-step', help='The step of the range.')

parser.add_argument('--re-download', help='Force re download of genomes.', action='store_true')
parser.add_argument('--refseq', help='Use RefSeq expanded search.', action='store_true')
parser.add_argument('--all-scores', help='Use all the tRNA-SeC, independent of score search.', action='store_false')

parser.add_argument('--quiet', help='Quiet printing.', action='store_false')
parser.add_argument('-QD', '--q-download', help='Quiet for download.', action='store_false')
parser.add_argument('-QF', '--q-fetch', help='Quiet for fetch.', action='store_false')
parser.add_argument('-QS', '--q-scan', help='Quiet for scan.', action='store_false')
parser.add_argument('-QDe', '--q-detected', help='Quiet for detected tRNAs-SeC.', action='store_false')
parser.add_argument('-QP', '--q-preprocess', help='Quiet for preprocess tRNAs-SeC.', action='store_false')
parser.add_argument('-QT', '--q-taxonomy', help='Quiet for taxonomy detection.', action='store_false')
parser.add_argument('-QTA', '--q-taxana', help='Quiet for taxonomy analysis.', action='store_false')
parser.add_argument('-QRRS', '--q-rRNAs', help='Quiet for rRNAsmallsubunit collection.', action='store_false')

parser.add_argument('-SD', '--suppress-download', help='Suppress download process.', action='store_false')
parser.add_argument('-SF', '--suppress-fetch', help='Suppress fetch process.', action='store_false')
parser.add_argument('-SS', '--suppress-scan', help='Suppress scan process.', action='store_false')
parser.add_argument('-SDe', '--suppress-detected', help='Suppress detected process.', action='store_false')
parser.add_argument('-SP', '--suppress-preprocess', help='Suppress preprocess.', action='store_false')
parser.add_argument('-SPh', '--suppress-taxonomy', help='Suppress taxonomy process.', action='store_false')
parser.add_argument('-SPA', '--suppress-taxana', help='Suppress taxonomy analysis process.', action='store_false')
parser.add_argument('-SRRS', '--suppress-rRNAs', help='Suppress rRNAsmallsubunit collection.', action='store_false')
parser.add_argument('-SAl', '--suppress-align', help='Suppress MAFFT alignment.', action='store_false')
parser.add_argument('-SM', '--suppress-metadata', help='Suppress metadata process.', action='store_false')

args=parser.parse_args()



genomesPath = args.genomes_path if args.genomes_path != None else 'GenomesArchaeaExpanded/'
speciesFile = args.species_file if args.species_file != None else 'species'
fetchFile = args.fetch_file if args.fetch_file != None else 'fetch'
readyFile = args.ready_file if args.ready_file != None else 'ready'
detectedFile = args.detected_file if args.detected_file != None else 'detected'
processedFile = args.processed_file if args.processed_file != None else 'processed'
alignFile = args.align_file if args.align_file != None else 'align'
taxonomyFile = args.taxonomy_file if args.taxonomy_file != None else 'taxonomy'
metadataFile = args.metadata_file if args.metadata_file != None else 'metadata'

read = args.read_species
save = args.save_species

taxonNames = eval(args.taxon_names) if args.taxon_names != None else __taxonNames
taxonLevel = args.taxon_level if args.taxon_level != None else 'all'
sequentialAnalysis = args.sequential

referenceRange = int(args.reference_range) if args.reference_range != None else None
rangeStep = int(args.range_step) if args.range_step != None else 1

reDownload = args.re_download
refseq = args.refseq
highestScore = args.all_scores

quiet = args.quiet
if not quiet:
    verboseDownload = 0
    verboseFetch = 0
    verboseScan = 0
    verboseDetected = 0
    verbosePreprocess = 0
    verboseTaxonomy = 0
    verboseTaxonAnalysis = 0
    verboseRRS = 0
else:
    verboseDownload = args.q_download
    verboseFetch = args.q_fetch
    verboseScan = args.q_scan
    verboseDetected = args.q_detected
    verbosePreprocess = args.q_preprocess
    verboseTaxonomy = args.q_taxonomy
    verboseTaxonAnalysis = args.q_taxana
    verboseRRS = args.q_rRNAs

suppressDownload = args.suppress_download
suppressFetch = args.suppress_fetch
suppressScan = args.suppress_scan
suppressDetected = args.suppress_detected
suppressPreprocess = args.suppress_preprocess
suppressTaxonomy = args.suppress_taxonomy
suppressTaxonAnalysis = args.suppress_taxana
suppressRRS = args.suppress_rRNAs
suppressAlign = args.suppress_align
suppressMetadata = args.suppress_metadata



suffix = ''
if referenceRange != None:
    suffix = f'_{referenceRange}*{rangeStep}'

suffixScore = ''
if not highestScore:
    suffixScore = '_all_scores'

files = {
    '__genomesPath': genomesPath,
    '__fetchFile': fetchFile + suffix,
    '__readyFile': readyFile + suffix,
    '__detectedFile': detectedFile + suffix,
    '__processedFile': processedFile + suffix + suffixScore,
    '__alignFile': alignFile + suffixScore,
    '__taxonomyFile': taxonomyFile,
    '__metadataFile': metadataFile + suffixScore,

    'suppressDownload': not suppressDownload,
    'suppressFetch': not suppressFetch,
    'suppressDetected': not suppressDetected,
    'suppressPreprocess': not suppressPreprocess,
    'suppressTaxonomy': not suppressTaxonomy,
    'suppressAlign': not suppressAlign
}
initiate(**files)

species = collectInfo(taxonNames, archaea=refseq, save=save, read=read, __speciesFile=speciesFile)

if suppressDownload: downloadGenomes(species, sizeLimit=20, referenceRange=referenceRange, rangeStep=rangeStep, reDownload=reDownload, verbose=verboseDownload)
if suppressFetch: downloadFetch(verbose=verboseFetch)
if suppressScan: trnaScanSE(verbose=verboseScan)
if suppressDetected: findDetectedSeC(verbose=verboseDetected)
if suppressTaxonomy: taxonCollection(verbose=verboseTaxonomy)
if suppressPreprocess: preprocessAndMetadata(highestScore=highestScore, verbose=verbosePreprocess, debug=0)
if suppressTaxonAnalysis: taxonAnalysis(taxonLevel, verbose=verboseTaxonAnalysis, sequential=sequentialAnalysis, debug=0)
# if suppressRRS: collectRRS(verbose=verboseRRS)

if ((referenceRange == None) and (suppressAlign)):
    alignMAFFT()