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
filelist=`find ../ucscdownfiltered_sort/ -iname '*_[1-2].bed'|xargs -l1 basename |grep '.bed' | cut -d '_' -f 1 | sort | uniq | grep Ag04450Uw`
parallel -j 20 --nice 19 workfunc {} ::: $filelist

# for f in $filelist;
# do
#     currfiles=(`find ../ucscdownfiltered/ -iname '*_[1-2].bed'| grep '.bed' | grep ${f}_`)  # parenthesis makes it an array
# #    echo `grep -c .bed $currfiles`
#     #echo $currfiles
#     # echo "${#currfiles[@]}"
#     if [ "${#currfiles[@]}" -eq 2 ]
#         then
# #            echo ${currfiles[@]}
# 	    python mergebedreps.py $f.bed ${currfiles[@]}
#     else
	
# #	echo $f
# #	echo ${currfiles[0]}
# 	sleep 1
# 	#cp ${currfiles[0]} $f.bed
#     fi
# done
