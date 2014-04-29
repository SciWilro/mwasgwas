
FNAME=$RANDOM
cut  -f 2- $1 > ${FNAME}
cut  -f 1 ${FNAME} > ${FNAME}1
paste ${FNAME}1 ${FNAME} > $1
rm ${FNAME}1 ${FNAME}
