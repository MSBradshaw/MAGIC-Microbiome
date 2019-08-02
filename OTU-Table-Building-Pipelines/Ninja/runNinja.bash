#!/bin/bash

rm *.fasta
while read FN
do
	file=$(echo $FN | cut -d'/' -f 3 | tr '\n' '.')
	folder=output/$file
	mkdir "$folder"
	suffix='fasta'
	suf='.fna'
	old=$file$suf
	newfile=$file$suffix
	sed -e 's/\(^>.*$\)/#\1#/' $FN | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > fastas/$newfile
	python ../../NINJA-OPS/bin/ninja.py -i fastas/$newfile -o $folder -b silva123_97
#	rm temp.fasta
	#python ../../NINJA-OPS/bin/ninja.py -i $newfile -o $file
done < <(find ../01.CleanData -maxdepth 2 -type f -name '*.fna')
#-printf "%f\n")
