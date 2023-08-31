cat *.abundance.txt | sed 's/#circ[[:digit:]]\+-/\t/' > realigned.csv
blastn -query BRAGUZ_M_07_Guzerat.srf.fa -subject BRAGUZ_M_07_Guzerat.srf.fa -outfmt 6 -best_hit_score_edge 0.1 | awk '$1<$2' | less
grep "BTSAT2\s" satellites/*out | awk '($8-$7)>1400' | grep -oP "\S+circ\S+" | sed 's/#circ[[:digit:]]\+-/\t/' | cut -f 2 | sort | uniq -c | sort -k1,1nr
