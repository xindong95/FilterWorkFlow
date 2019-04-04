for dir in ./GSM*
do
    echo $dir
    wc -l ./${dir}/${dir}_peaks.narrowPeak >> ./peaks_summary.txt
done
echo ========================================= >> ./peaks_summary.txt
for file in ./10FoldChange/*
do
    wc -l ${file} >> ./peaks_summary.txt
done
echo ========================================= >> ./peaks_summary.txt
for file in ./20FoldChange/*
do
    wc -l ${file} >> ./peaks_summary.txt
done
