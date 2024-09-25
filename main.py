from datasets import collectInfo, downloadGenomes, downloadFetch, trnaScanSE, findDetectedSeC, preprocessSeC, alignMAFFT
from datasets import initiate, pretty
import os, argparse



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
parser.add_argument('--range-step', help='The step of the range')
args=parser.parse_args()

genomesPath = args.genomes_path if args.genomes_path != None else 'Genomes/'
fetchFile = args.fetch_file if args.fetch_file != None else 'fetch'
readyFile = args.ready_file if args.ready_file != None else 'ready'
detectedFile = args.detected_file if args.detected_file != None else 'detected'
processedFile = args.processed_file if args.processed_file != None else 'processed'
alignFile = args.align_file if args.align_file != None else 'align'
taxons = eval(args.taxon) if args.taxon != None else __taxons
referenceRange = int(args.reference_range) if args.taxon != None else None
rangeStep = int(args.range_step) if args.range_step != None else 1



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

species = collectInfo(taxons)
downloadGenomes(species, sizeLimit=20, referenceRange=referenceRange, rangeStep=rangeStep)
downloadFetch()
trnaScanSE()
findDetectedSeC()
preprocessSeC()

if referenceRange == None:
    alignMAFFT()