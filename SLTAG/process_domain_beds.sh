#/bin/env bash

# grep domain lists
grep -f domain_list_LTR_only.txt cdd.bed | sort -k1,1 -k2,2n > LINE.0
grep -f domain_list_LINE_only.txt cdd.bed | sort -k1,1 -k2,2n > LTR.0
grep -f domain_list_LTR_or_LINE.txt cdd.bed | sort -k1,1 -k2,2n > LTR_or_LINE.0
grep -f domain_list_not_LINE_or_LTR.txt cdd.bed | sort -k1,1 -k2,2n > other.0

# remove any overlaps from ‘other’
bedtools subtract -a other.0 -b LTR.0 > other.1
bedtools subtract -a other.1 -b LINE.0 > other.2
bedtools subtract -a other.2 -b LTR_or_LINE.0 > other.final

# remove any LINE-only or LTR-only from LTR-or-LINE
bedtools subtract -a LTR_or_LINE.0 -b LINE.0 > LTR_or_LINE.1
bedtools subtract -a LTR_or_LINE.1 -b LTR.0 > LTR_or_LINE.2

# if there are any segments in both LTR-only and LINE-only, remove from both and add to LTR-or-LINE
bedtools intersect -a LTR.0 -b LINE.0  >LTR_intersect_LINE.0
bedtools subtract -a LTR.0 -b LTR_intersect_LINE.0 > LTR.final
bedtools subtract -a LINE.0 -b LTR_intersect_LINE.0 > LINE.final
bedtools merge -a LTR_or_LINE.2 -b LTR_intersect_LINE.0 > LTR_or_LINE.final

# name the features in each file and concatenate them
sed -i 's/$/\t0\t+\tLTR/‘ LTR.final
sed -i 's/$/\t0\t+\tLINE/‘ LINE.final
sed -i 's/$/\t0\t+\tLTR|LINE/‘ LTR_or_LINE.final
sed -i 's/$/\t0\t+\tother/‘ other.final
cat *.final | sort -k1,1 -k2,2n >TE_domains.bed

echo Number overlapping (should be 0):
bedtools merge -n TE_domains.bed

