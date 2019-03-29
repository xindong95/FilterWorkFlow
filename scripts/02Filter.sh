for file in ./sorted_bam/*.bam
do
    filename=`basename $file`
    echo $filename
    mkdir ./${filename}
    samtools view -q 1 -@ 8 $file | cut -f9 > ./${filename}/fragment.txt
    echo 'plot distribution'
    Rscript ./script/density.r ./${filename}/fragment.txt ./${filename}/${filename}_distribution.png
    # echo 'start filter'
    # mkdir ./${filename}/bam_file/
    # samtools view -h $file | awk '($9 <= 50 && $9 >= -50) || $1 ~ /^@/' | samtools view -bS - > ./${filename}/bam_file/${filename}.sub50.bam
    # samtools view -h $file | awk '($9 <= 100 && $9 >= -100) || $1 ~ /^@/' | samtools view -bS - > ./${filename}/bam_file/${filename}.sub100.bam
    # samtools view -h $file | awk '($9 <= 150 && $9 >= -150) || $1 ~ /^@/' | samtools view -bS - > ./${filename}/bam_file/${filename}.sub150.bam
    # samtools view -h $file | awk '($9 <= 200 && $9 >= -200) || $1 ~ /^@/' | samtools view -bS - > ./${filename}/bam_file/${filename}.sub200.bam
    # samtools view -h $file | awk '($9 <= 250 && $9 >= -250) || $1 ~ /^@/' | samtools view -bS - > ./${filename}/bam_file/${filename}.sub250.bam
    # samtools view -h $file | awk '($9 <= 300 && $9 >= -300) || $1 ~ /^@/' | samtools view -bS - > ./${filename}/bam_file/${filename}.sub300.bam
    # samtools view -h $file | awk '($9 <= 350 && $9 >= -350) || $1 ~ /^@/' | samtools view -bS - > ./${filename}/bam_file/${filename}.sub350.bam
    # samtools view -h $file | awk '($9 <= 400 && $9 >= -400) || $1 ~ /^@/' | samtools view -bS - > ./${filename}/bam_file/${filename}.sub400.bam
    # samtools view -h $file | awk '($9 <= 450 && $9 >= -450) || $1 ~ /^@/' | samtools view -bS - > ./${filename}/bam_file/${filename}.sub450.bam
    # samtools view -h $file | awk '($9 <= 500 && $9 >= -500) || $1 ~ /^@/' | samtools view -bS - > ./${filename}/bam_file/${filename}.sub500.bam
done

