#!/bin/sh

if [ $# != 1 ]
then
    echo "USAGE: $0 matched_set.gtf"
    exit 255;
fi

echo "PASA-transcript count:"
grep transcript $1 | grep -c start_
#grep cufflinks-transcript $1 | grep -c start_
echo "Exonerate count:"
grep Exonerate $1 | grep -c start_
echo "TransposonPSI count:"
grep TransposonPSI $1 | grep -c start_
