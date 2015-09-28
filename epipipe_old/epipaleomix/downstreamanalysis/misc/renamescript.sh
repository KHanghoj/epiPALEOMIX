find . -type f | while read FILE
do
    # modify line below to do what you need, then remove leading "echo"
    #echo mv $FILE $(echo $FILE | sed 's/LBK//g')
    mv $FILE $(echo $FILE | sed 's/Kostenki/k14/g')
done
