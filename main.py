from datasets import collectInfo, downloadGenomes, downloadFetch, trnaScanSE, findDetectedSeC, preprocessSeC
from datasets import initiate, pretty
import os, argparse



os.system('clear')



parser=argparse.ArgumentParser()
parser.add_argument("--genomes-path", help="Path for the genomes to be stored.")
parser.add_argument("--fetch-file", help="The name for the default fetch file information.")
parser.add_argument("--ready-file", help="The name for the default ready file information.")
parser.add_argument("--detected-file", help="The name for the default detected file information.")
parser.add_argument("--processed-file", help="The name for the default processed file information.")
args=parser.parse_args()



genomesPath = args.genomes_path if args.genomes_path != None else 'Genomes/'
fetchFile = args.fetch_file if args.fetch_file != None else 'fetch'
readyFile = args.ready_file if args.ready_file != None else 'ready'
detectedFile = args.detected_file if args.detected_file != None else 'detected'
processedFile = args.processed_file if args.processed_file != None else 'processed'



taxons = {
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

initiate(__genomesPath=genomesPath, __fetchFile=fetchFile, __readyFile=readyFile, __detectedFile=detectedFile, __processedFile=processedFile)

species = collectInfo(taxons)
downloadGenomes(species, sizeLimit=20)
downloadFetch()
trnaScanSE()
findDetectedSeC()
preprocessSeC()