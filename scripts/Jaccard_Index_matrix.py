import os
import subprocess
import sys
from optparse import OptionParser

def getJaccard(peaks_a,peaks_b):
    peaks_a = peaks_a #'GSM2644567.raw/GSM2644567.raw_peaks.bed'
    peaks_b = peaks_b #'GSM2644567.sub150/GSM2644567.sub150_peaks.bed'
    cmd = ['bedtools','jaccard','-a', peaks_a,'-b', peaks_b]
    jaccard = subprocess.check_output(cmd).decode('utf8').split('\n')[1].split('\t')[2]
    # jaccard = float(jaccard_str)
    return jaccard

if __name__ == '__main__':
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-d", "--dir", help="dir")
    optparser.add_option("-o", "--output", help="output matrix files")
    optparser.add_option("-r", "--rscript", help="rscript")
    optparser.add_option("-p", "--plot", help="plot")
    (options, args) = optparser.parse_args(sys.argv)

    directory = os.path.abspath(options.dir)
    folders = []
    for i in os.listdir(directory):
        if i.startswith('GSM'):
            folders.append(i)
        else:
            continue

    folders.sort()
    header = ['samples']
    header.extend(folders)

    matrix = [header]

    for i in folders:
        peaks_a = ''.join([directory,'/',i,'/',i,'_peaks.bed'])
        row = [i]
        for j in folders:
            peaks_b = ''.join([directory,'/',j,'/',j,'_peaks.bed'])
            row.append(getJaccard(peaks_a,peaks_b))
        matrix.append(row)
    table = []
    for i in matrix:
        table.append('\t'.join(i))
    with open(options.output, 'w') as output:
        output.write('\n'.join(table))

    plot_cmd = ['Rscript', options.rscript, options.output, options.plot]
    subprocess.check_output(plot_cmd)



