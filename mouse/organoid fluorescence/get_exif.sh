#!/bin/bash
# find all the tif files (recursievely) 
# use exiftool to extract metadata 
# and save to a file removing .tif and appending _meta.txt

for f in $(find . -name '*.tif')
do
    echo "$f"
    # echo "${f%.tif}""_meta.txt"
    exiftool $f > "${f%.tif}""_meta.txt"
done
