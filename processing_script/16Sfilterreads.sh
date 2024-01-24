for i in *.fq
do
/WorkingDogs/FilterReads/ubuntu-18.04/FilterReads -r 16S -t 8 -full -qt 30 +fz ~/cohort/CAMI/WorkingDogs/BuildFilter/ubuntu-18.04/Silva-RefSeq_5-18_16S_NC_20.mer 50% ~/cohort/CAMI/16S/$i
done
