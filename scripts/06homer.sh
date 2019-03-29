source activate chips
for file in ./sorted_bam/*.bam
do
    filename=`basename $file`
    echo $filename
    findMotifsGenome.pl ${filename}/peaks/none_filter/${filename}_None_filter_summits.bed hg38 ${filename}/motif/none_filter/ -size 600 -p 8 -mask -seqlogo -preparsedDir ${filename}/motif/none_filter/
    findMotifsGenome.pl ${filename}/peaks/sub50/${filename}_sub50_summits.bed hg38 ${filename}/motif/sub50/ -size 600 -p 8 -mask -seqlogo -preparsedDir ${filename}/motif/sub50/
    findMotifsGenome.pl ${filename}/peaks/sub100/${filename}_sub100_summits.bed hg38 ${filename}/motif/sub100/ -size 600 -p 8 -mask -seqlogo -preparsedDir ${filename}/motif/sub100/
    findMotifsGenome.pl ${filename}/peaks/sub150/${filename}_sub150_summits.bed hg38 ${filename}/motif/sub150/ -size 600 -p 8 -mask -seqlogo -preparsedDir ${filename}/motif/sub150/
    findMotifsGenome.pl ${filename}/peaks/sub200/${filename}_sub200_summits.bed hg38 ${filename}/motif/sub200/ -size 600 -p 8 -mask -seqlogo -preparsedDir ${filename}/motif/sub200/
done
conda deactivate