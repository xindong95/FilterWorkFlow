for dir in /mnt/Storage/home/dongxin/Storage2/Workflow_ATAC_Filter/analysis/GSM*
do
	echo $dir
    python ~/Projects/FilterWorkFlow/scripts/Jaccard_Index_matrix.py -d $dir/peaks/ -o $dir/peaks/jaccard_matrix.txt -r ~/Projects/FilterWorkFlow/scripts/Jaccard_Index_heatmap.r -p $dir/peaks/jaccard.png
done