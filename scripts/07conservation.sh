for file in ./sorted_bam/*.bam
do
    filename=`basename $file`
    echo ${filename}
    head -n 5000 ./${filename}/peaks/none_filter/${filename}_None_filter_summits.bed > ./${filename}/peaks/none_filter/${filename}_None_filter_sorted_5000_summits.bed
    head -n 5000 ./${filename}/peaks/sub50/${filename}_sub50_summits.bed > ./${filename}/peaks/sub50/${filename}_sub50_sorted_5000_summits.bed
    head -n 5000 ./${filename}/peaks/sub100/${filename}_sub100_summits.bed > ./${filename}/peaks/sub100/${filename}_sub100_sorted_5000_summits.bed
    head -n 5000 ./${filename}/peaks/sub150/${filename}_sub150_summits.bed > ./${filename}/peaks/sub150/${filename}_sub150_sorted_5000_summits.bed
    head -n 5000 ./${filename}/peaks/sub200/${filename}_sub200_summits.bed > ./${filename}/peaks/sub200/${filename}_sub200_sorted_5000_summits.bed
    mkdir ./${filename}/conserv/
    mkdir ./${filename}/conserv/none_filter/
    mkdir ./${filename}/conserv/sub50/
    mkdir ./${filename}/conserv/sub100/
    mkdir ./${filename}/conserv/sub150/
    mkdir ./${filename}/conserv/sub200/
    python script/conservation_plot.py -t Conservation_at_summits -d ~/Files/ref_files/mm10/conservation/vertebrate -o ./${filename}/conserv/none_filter/plot -l Peak_summits ./${filename}/peaks/none_filter/${filename}_None_filter_sorted_5000_summits.bed -w 4000 > ./${filename}/conserv/none_filter/conserv.txt
    python script/conservation_plot.py -t Conservation_at_summits -d ~/Files/ref_files/mm10/conservation/vertebrate -o ./${filename}/conserv/sub50/plot -l Peak_summits ./${filename}/peaks/sub50/${filename}_sub50_sorted_5000_summits.bed -w 4000 > ./${filename}/conserv/sub50/conserv.txt
    python script/conservation_plot.py -t Conservation_at_summits -d ~/Files/ref_files/mm10/conservation/vertebrate -o ./${filename}/conserv/sub100/plot -l Peak_summits ./${filename}/peaks/sub100/${filename}_sub100_sorted_5000_summits.bed -w 4000 > ./${filename}/conserv/sub100/conserv.txt
    python script/conservation_plot.py -t Conservation_at_summits -d ~/Files/ref_files/mm10/conservation/vertebrate -o ./${filename}/conserv/sub150/plot -l Peak_summits ./${filename}/peaks/sub150/${filename}_sub150_sorted_5000_summits.bed -w 4000 > ./${filename}/conserv/sub150/conserv.txt
    python script/conservation_plot.py -t Conservation_at_summits -d ~/Files/ref_files/mm10/conservation/vertebrate -o ./${filename}/conserv/sub200/plot -l Peak_summits ./${filename}/peaks/sub200/${filename}_sub200_sorted_5000_summits.bed -w 4000 > ./${filename}/conserv/sub200/conserv.txt
done