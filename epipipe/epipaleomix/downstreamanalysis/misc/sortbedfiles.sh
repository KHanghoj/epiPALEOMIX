function dosort(){
    echo $1
    f=$1
    sort -k 1,1n -k2,2n -k 3,3n $f |sed '/^X/ d' | sed '/^Y/ d'| sed '/^M/ d' > ${f}.tmp
    cp ${f}.tmp $f
    awk '{print "chr" $0;}' ${f}.tmp > ${f}.tmp.tmp && mv ${f}.tmp.tmp ${f}.tmp
    rename 's/_wochr.bed.tmp/_wchr.bed/' ${f}.tmp
    
}

files=`find /home/krishang/data/bedfiles -iname "*wochr.bed" | grep -v oldstuff | grep 2000`
export -f dosort
parallel -j 5 --nice 19 dosort {} ::: $files



