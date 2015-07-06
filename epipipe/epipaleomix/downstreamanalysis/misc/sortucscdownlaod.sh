#!/bin/bash
function workfuncsort(){
    f=$1
    #echo $f
    #echo ${f##*/}
    sort -k 1,1n -k 2,2n $f > ${f##*/}
}
export -f workfuncsort
filelist=`find ../ucscdownfiltered/ -iname '*_[1-2].bed' | grep .bed`
parallel -j 20 --nice 19 workfuncsort {} ::: $filelist

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
