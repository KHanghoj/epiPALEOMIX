#!/bin/bash
function workfunc(){
    f=$1
    currfiles=(`find ../ucscdownfiltered_sort/ -iname '*_[1-2].bed'| grep '.bed' | grep ${f}_`)  # parenthesis makes it an array
    if [ "${#currfiles[@]}" -eq 2 ]; then
	python mergebedreps.py $f.bed ${currfiles[@]}
    else  # in case only one, then just copied.
	echo "Just a single file: ${currfiles[@]}"
	cp ${currfiles[0]} $f.bed
    fi
}
export -f workfunc
filelist=`find ../ucscdownfiltered_sort/ -iname '*_[1-2].bed'|xargs -l1 basename |grep '.bed' | cut -d '_' -f 1 | sort | uniq `
parallel -j 20 --nice 19 workfunc {} ::: $filelist
