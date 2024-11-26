#!/bin/bash

source ~/.bashrc
source activate base

TRNA="treefile/align.fasta.treefile"
RRNA="treefile/alignRSSU.fasta.treefile"

echo ""
echo "tqDist methods"
echo "Quartet distance:"
echo -e "#leaves\t#quartets\tdist\t\tnormDist\t#resQuartetsA\t<--norm\t\t#unrQua\t<--norm"
quartet_dist -v "$RRNA" "$TRNA"

echo ""
echo "Pairs quartet distance:"
echo -e "#leaves\t#quartets\tdist\t\tnormDist\t#resQuartetsA\t<--norm\t\t#unrQua\t<--norm"
pairs_quartet_dist -v "$RRNA" "$TRNA"

echo ""
echo "Triplet distance:"
echo -e "#leaves\t#triplets\tdist\t\tnormDist\t#resTripletsA\t<--norm\t\t#unrTri\t<--norm"
triplet_dist -v "$RRNA" "$TRNA"

echo ""
echo "Pair triplet distance:"
echo -e "#leaves\t#triplets\tdist\t\tnormDist\t#resTripletsA\t<--norm\t\t#unrTri\t<--norm"
pairs_triplet_dist -v "$RRNA" "$TRNA"

echo ""
echo "ETE3 Compare method (Robinson–Foulds):"
conda activate ete3
ete3 compare -r "$RRNA" -t "$TRNA" --unrooted
conda deactivate

echo ""
echo "Python Kendall–Colijn distance (lam = .85):"
python3 -c "from treefile import KC_dist;import os;res = KC_dist.KC_dist(os.popen('cat $RRNA').read().strip(), os.popen('cat $TRNA').read().strip(), lam=.85);print(res)"
echo "Python Kendall–Colijn distance (lam = .15):"
python3 -c "from treefile import KC_dist;import os;res = KC_dist.KC_dist(os.popen('cat $RRNA').read().strip(), os.popen('cat $TRNA').read().strip(), lam=.15);print(res)"
