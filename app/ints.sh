#!/bin/bash
gmsfile=$1
exotimofile=$2
sed -n '/2 ELECTRON/,/TOTAL NUMBER/p' $gmsfile | awk '$1 ~ /[0-9]/ && $2 !~ /ELECTRON/' > tmp1.$$
awk '{print $1,$2,$3,$4,$6; print $7,$8,$9,$10,$12;}' tmp1.$$ > tmp2.$$
sort < tmp2.$$ > gms.ints
sort < $exotimofile > exotimo.ints
diff gms.ints exotimo.ints > ints.diff
rm *.$$
