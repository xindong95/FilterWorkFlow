source activate chips
mkdir bam
mkdir sorted_bam
for file in ./fastq/*_R1.fastq
do
	echo ${file}
	filename=${file##*/}
	GSM=${filename%_*}
    bwa mem -t 8 ~/Files/ref_files/hg38/bwa_indices/hg38/hg38.fa ${file} ${file/R1/R2} | samtools view -Sb - > ./bam/${GSM}.bam
    sambamba sort ./bam/${GSM}.bam -o ./sorted_bam/${GSM} -t 8
done
conda deactivate