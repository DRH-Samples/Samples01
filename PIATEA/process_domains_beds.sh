#/bin/env bash

echo Grep domain lists...
grep -f domain_list_LTR_only.txt cdd.bed | sort -k1,1 -k2,2n > LTR.0
grep -f domain_list_LINE_only.txt cdd.bed | sort -k1,1 -k2,2n > LINE.0
grep -f domain_list_LTR_or_LINE.txt cdd.bed | sort -k1,1 -k2,2n > LTR_or_LINE.0
grep -f domain_list_not_LINE_or_LTR.txt cdd.bed | sort -k1,1 -k2,2n > other.0

echo Remove overlaps from 'other'...
bedtools subtract -a other.0 -b LTR.0 > other.1
bedtools subtract -a other.1 -b LINE.0 > other.2
bedtools subtract -a other.2 -b LTR_or_LINE.0 > other.3

echo Remove any LINE-only or LTR-only from LTR-or-LINE...
bedtools subtract -a LTR_or_LINE.0 -b LINE.0 > LTR_or_LINE.1
bedtools subtract -a LTR_or_LINE.1 -b LTR.0 > LTR_or_LINE.2

echo If there are any segments in both LTR-only and LINE-only, remove from both and add to LTR-or-LINE...
bedtools intersect -a LTR.0 -b LINE.0  >LTR_intersect_LINE.0
if [ -s LTR_intersect_LINE.0 ]
then
	bedtools subtract -a LTR.0 -b LTR_intersect_LINE.0 > LTR.1
	bedtools subtract -a LINE.0 -b LTR_intersect_LINE.0 > LINE.1
	bedtools merge -a LTR_or_LINE.2 -b LTR_intersect_LINE.0 > LTR_or_LINE.3
else
	cp LTR.0 LTR.1
	cp LINE.0 LINE.1
	cp LTR_or_LINE.2 LTR_or_LINE.3
fi

echo Flatten, name, and concatenate...
bedtools sort -i LTR.1 | bedtools merge > LTR.final
bedtools sort -i LINE.1 | bedtools merge > LINE.final
bedtools sort -i LTR_or_LINE.3 | bedtools merge > LTR_or_LINE.final
bedtools sort -i other.3 | bedtools merge > other.final
sed -i 's/$/\tLTR/' LTR.final
sed -i 's/$/\tLINE/' LINE.final
sed -i 's/$/\tLTR|LINE/' LTR_or_LINE.final
sed -i 's/$/\tother/' other.final
cat *.final | bedtools sort >TE_domains.bed

echo Number of records that overlap by at least 1 bp \(should be 0\)...
bedtools merge -n -d -1 -i TE_domains.bed | grep -v -e '1$' | wc -l

