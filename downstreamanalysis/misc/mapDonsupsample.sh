#!/bin/bash

function docalculations () {
    file=$1
    x=`samtools view -c $file`
    prop=`python -c "print 300000.0/$x"`
    samtools view -s  $prop -h -b $1 > subsample_${file##*/}
    samtools index subsample_${file##*/}
    samtools view -h -b  subsample_${file##*/} '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' > subsample_keep_${file##*/}
    rm subsample_${file##*/}
    rm subsample_${file##*/}.bai
    samtools view -H subsample_keep_${file##*/}| sed '/^@SQ\tSN:GL/ d' | samtools reheader - subsample_keep_${file##*/} > subsample_final_${file##*/}
    rm subsample_keep_${file##*/}
}
export -f docalculations
files=`find ../ -iname '*bam' | grep -v 1000g`
echo $files
parallel -j 20 --nice 19 docalculations {} ::: $files


for file in `find . -iname '*bam' | grep _final_`; do echo "nice -n 19 mapDamage -n 100000 -i ${file} -r /home/krishang/data/reference_human/hs.build37.1.fa –burn=25000 –iter=100000 --no-stats"; done | parallel -j 20
