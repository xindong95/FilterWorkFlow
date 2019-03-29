source activate chips
for file in ./sorted_bam/*.bam
do
    filename=`basename $file`
    dir=`pwd`
    echo $filename
    cd ${filename}/peaks/
    wc -l ./none_filter/${filename}_None_filter_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./sub50/${filename}_sub50_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./sub100/${filename}_sub100_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./sub150/${filename}_sub150_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./sub200/${filename}_sub200_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./sub250/${filename}_sub250_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./sub300/${filename}_sub300_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./sub350/${filename}_sub350_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./sub400/${filename}_sub400_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./sub450/${filename}_sub450_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./sub500/${filename}_sub500_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    echo =================================================== >> ${dir}/peaks_summary.txt

    mkdir 10FoldChange
    awk '$7 > 10' ./none_filter/${filename}_None_filter_peaks.narrowPeak > ./10FoldChange/${filename}_None_filter_peaks.narrowPeak
    awk '$7 > 10' ./sub50/${filename}_sub50_peaks.narrowPeak > ./10FoldChange/${filename}_sub50_peaks.narrowPeak
    awk '$7 > 10' ./sub100/${filename}_sub100_peaks.narrowPeak > ./10FoldChange/${filename}_sub100_peaks.narrowPeak
    awk '$7 > 10' ./sub150/${filename}_sub150_peaks.narrowPeak > ./10FoldChange/${filename}_sub150_peaks.narrowPeak
    awk '$7 > 10' ./sub200/${filename}_sub200_peaks.narrowPeak > ./10FoldChange/${filename}_sub200_peaks.narrowPeak
    # awk '$7 > 10' ./sub250/${filename}_sub250_peaks.narrowPeak > ./10FoldChange/${filename}_sub250_peaks.narrowPeak
    # awk '$7 > 10' ./sub300/${filename}_sub300_peaks.narrowPeak > ./10FoldChange/${filename}_sub300_peaks.narrowPeak
    # awk '$7 > 10' ./sub350/${filename}_sub350_peaks.narrowPeak > ./10FoldChange/${filename}_sub350_peaks.narrowPeak
    # awk '$7 > 10' ./sub400/${filename}_sub400_peaks.narrowPeak > ./10FoldChange/${filename}_sub400_peaks.narrowPeak
    # awk '$7 > 10' ./sub450/${filename}_sub450_peaks.narrowPeak > ./10FoldChange/${filename}_sub450_peaks.narrowPeak
    # awk '$7 > 10' ./sub500/${filename}_sub500_peaks.narrowPeak > ./10FoldChange/${filename}_sub500_peaks.narrowPeak

    wc -l ./10FoldChange/${filename}_None_filter_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./10FoldChange/${filename}_sub50_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./10FoldChange/${filename}_sub100_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./10FoldChange/${filename}_sub150_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./10FoldChange/${filename}_sub200_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./10FoldChange/${filename}_sub250_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./10FoldChange/${filename}_sub300_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./10FoldChange/${filename}_sub350_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./10FoldChange/${filename}_sub400_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./10FoldChange/${filename}_sub450_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./10FoldChange/${filename}_sub500_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    echo =================================================== >> ${dir}/peaks_summary.txt

    mkdir 20FoldChange
    awk '$7 > 20' ./none_filter/${filename}_None_filter_peaks.narrowPeak > ./20FoldChange/${filename}_None_filter_peaks.narrowPeak
    awk '$7 > 20' ./sub50/${filename}_sub50_peaks.narrowPeak > ./20FoldChange/${filename}_sub50_peaks.narrowPeak
    awk '$7 > 20' ./sub100/${filename}_sub100_peaks.narrowPeak > ./20FoldChange/${filename}_sub100_peaks.narrowPeak
    awk '$7 > 20' ./sub150/${filename}_sub150_peaks.narrowPeak > ./20FoldChange/${filename}_sub150_peaks.narrowPeak
    awk '$7 > 20' ./sub200/${filename}_sub200_peaks.narrowPeak > ./20FoldChange/${filename}_sub200_peaks.narrowPeak
    # awk '$7 > 20' ./sub250/${filename}_sub250_peaks.narrowPeak > ./20FoldChange/${filename}_sub250_peaks.narrowPeak
    # awk '$7 > 20' ./sub300/${filename}_sub300_peaks.narrowPeak > ./20FoldChange/${filename}_sub300_peaks.narrowPeak
    # awk '$7 > 20' ./sub350/${filename}_sub350_peaks.narrowPeak > ./20FoldChange/${filename}_sub350_peaks.narrowPeak
    # awk '$7 > 20' ./sub400/${filename}_sub400_peaks.narrowPeak > ./20FoldChange/${filename}_sub400_peaks.narrowPeak
    # awk '$7 > 20' ./sub450/${filename}_sub450_peaks.narrowPeak > ./20FoldChange/${filename}_sub450_peaks.narrowPeak
    # awk '$7 > 20' ./sub500/${filename}_sub500_peaks.narrowPeak > ./20FoldChange/${filename}_sub500_peaks.narrowPeak

    wc -l ./20FoldChange/${filename}_None_filter_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./20FoldChange/${filename}_sub50_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./20FoldChange/${filename}_sub100_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./20FoldChange/${filename}_sub150_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    wc -l ./20FoldChange/${filename}_sub200_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./20FoldChange/${filename}_sub250_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./20FoldChange/${filename}_sub300_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./20FoldChange/${filename}_sub350_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./20FoldChange/${filename}_sub400_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./20FoldChange/${filename}_sub450_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    # wc -l ./20FoldChange/${filename}_sub500_peaks.narrowPeak >> ${dir}/peaks_summary.txt
    echo =================================================== >> ${dir}/peaks_summary.txt
    echo >> ${dir}/peaks_summary.txt
    echo =================================================== >> ${dir}/peaks_summary.txt
    cd ${dir}/
done
conda deactivate
