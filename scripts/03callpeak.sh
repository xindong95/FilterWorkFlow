source activate chips_py2
# for file in ./sorted_bam/*.bam
# do
#     filename=`basename $file`
#     echo $filename
echo 'start callpeak'
python2 ~/Applications/miniconda3/envs/chips_py2/bin/macs2 callpeak -B -q 0.01 --keep-dup 1 -g hs -f BAMPE --extsize 146 --nomodel -t ./sorted_bam/GSM1071280.sorted.bam -c ./sorted_bam/GSM1071276.sorted.bam --outdir ./GSM1071280.sorted.bam/peaks/none_filter/ -n GSM1071280.sorted.bam_None_filter
python2 ~/Applications/miniconda3/envs/chips_py2/bin/macs2 callpeak -B -q 0.01 --keep-dup 1 -g hs -f BAMPE --extsize 146 --nomodel -t ./GSM1071280.sorted.bam/bam_file/GSM1071280.sorted.bam.sub50.bam -c ./sorted_bam/GSM1071276.sorted.bam --outdir ./GSM1071280.sorted.bam/peaks/sub50/ -n GSM1071280.sorted.bam_sub50
python2 ~/Applications/miniconda3/envs/chips_py2/bin/macs2 callpeak -B -q 0.01 --keep-dup 1 -g hs -f BAMPE --extsize 146 --nomodel -t ./GSM1071280.sorted.bam/bam_file/GSM1071280.sorted.bam.sub100.bam -c ./sorted_bam/GSM1071276.sorted.bam --outdir ./GSM1071280.sorted.bam/peaks/sub100/ -n GSM1071280.sorted.bam_sub100
python2 ~/Applications/miniconda3/envs/chips_py2/bin/macs2 callpeak -B -q 0.01 --keep-dup 1 -g hs -f BAMPE --extsize 146 --nomodel -t ./GSM1071280.sorted.bam/bam_file/GSM1071280.sorted.bam.sub150.bam -c ./sorted_bam/GSM1071276.sorted.bam --outdir ./GSM1071280.sorted.bam/peaks/sub150/ -n GSM1071280.sorted.bam_sub150
python2 ~/Applications/miniconda3/envs/chips_py2/bin/macs2 callpeak -B -q 0.01 --keep-dup 1 -g hs -f BAMPE --extsize 146 --nomodel -t ./GSM1071280.sorted.bam/bam_file/GSM1071280.sorted.bam.sub200.bam -c ./sorted_bam/GSM1071276.sorted.bam --outdir ./GSM1071280.sorted.bam/peaks/sub200/ -n GSM1071280.sorted.bam_sub200

conda deactivate