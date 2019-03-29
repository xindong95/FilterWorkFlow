source activate chips
for file in ./sorted_bam/*.bam
do
	filename=`basename $file`
    echo $filename
	echo $filename >> ./intersect_peaks.txt
	wc -l ${filename}/peaks/sub150/${filename}_sub150_peaks.narrowPeak >> ./intersect_peaks.txt
	wc -l ${filename}/peaks/none_filter/${filename}_None_filter_peaks.narrowPeak >> ./intersect_peaks.txt
    bedtools intersect -wa -a ${filename}/peaks/sub150/${filename}_sub150_peaks.narrowPeak -b ${filename}/peaks/none_filter/${filename}_None_filter_peaks.narrowPeak > ${filename}/peaks/${filename}_intersect_peaks.narrowPeak
    bedtools intersect -v -a ${filename}/peaks/sub150/${filename}_sub150_peaks.narrowPeak -b ${filename}/peaks/none_filter/${filename}_None_filter_peaks.narrowPeak > ${filename}/peaks/${filename}_additional_peaks.narrowPeak
    wc -l ${filename}/peaks/${filename}_intersect_peaks.narrowPeak >> ./intersect_peaks.txt
    wc -l ${filename}/peaks/${filename}_additional_peaks.narrowPeak >> ./intersect_peaks.txt
    bedtools intersect -wa -a ${filename}/peaks/sub150/${filename}_sub150_peaks.narrowPeak -b /mnt/Storage/home/dongxin/Files/ref_files/hg38/regions/DHS_hg38.bed > ${filename}/peaks/${filename}_sub150_DHS_peaks.narrowPeak
    bedtools intersect -wa -a ${filename}/peaks/none_filter/${filename}_None_filter_peaks.narrowPeak -b /mnt/Storage/home/dongxin/Files/ref_files/hg38/regions/DHS_hg38.bed > ${filename}/peaks/${filename}_none_filter_DHS_peaks.narrowPeak
    bedtools intersect -wa -a ${filename}/peaks/${filename}_additional_peaks.narrowPeak -b /mnt/Storage/home/dongxin/Files/ref_files/hg38/regions/DHS_hg38.bed > ${filename}/peaks/${filename}_additional_DHS_peaks.narrowPeak
    echo ===DHS=== >> ./intersect_peaks.txt
    wc -l ${filename}/peaks/${filename}_sub150_DHS_peaks.narrowPeak >> ./intersect_peaks.txt
    wc -l ${filename}/peaks/${filename}_none_filter_DHS_peaks.narrowPeak >> ./intersect_peaks.txt
    wc -l ${filename}/peaks/${filename}_additional_DHS_peaks.narrowPeak >> ./intersect_peaks.txt
    echo ================================= >> ./intersect_peaks.txt
done
conda deactivate