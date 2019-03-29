#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import yaml
import math
from string import Template

# def getRuns(config):
#     """parse metasheet for Run groupings"""
#     ret = {}

#     #LEN: Weird, but using pandas to handle the comments in the file
#     #KEY: need skipinitialspace to make it fault tolerant to spaces!
#     metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#', skipinitialspace=True)
#     f = metadata.to_csv().split() #make it resemble an actual file with lines
#     #SKIP the hdr
#     for l in f[1:]:
#         tmp = l.strip().split(",")
#         #print(tmp)
#         ret[tmp[0]] = tmp[1:]

#     #print(ret)
#     config['runs'] = ret
#     return config

def addPy2Paths_Config(config):
    """ADDS the python2 paths to config"""
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    config["python2_pythonpath"] = os.path.join(conda_root, 'envs', 'chips_py2', 'lib', 'python2.7', 'site-packages')
    
    if not "python2" in config or not config["python2"]:
        config["python2"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'python2.7')

    if not "mdseqpos_path" in config or not config["mdseqpos_path"]:
        config["mdseqpos_path"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'MDSeqPos.py')

    if not "macs2_path" in config or not config["macs2_path"]:
        config["macs2_path"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'macs2')

def loadRef(config):
    """Adds the static reference paths found in config['ref']
    NOTE: if the elm is already defined, then we DO NOT clobber the value
    """
    f = open(config['ref'])
    ref_info = yaml.safe_load(f)
    f.close()
    #print(ref_info[config['assembly']])
    for (k,v) in ref_info[config['assembly']].items():
        #NO CLOBBERING what is user-defined!
        if k not in config:
            config[k] = v

#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config 
# config = getRuns(config)
addPy2Paths_Config(config)

#NOW load ref.yaml - SIDE-EFFECT: loadRef CHANGES config
loadRef(config)
#-----------------------------------------

# #------------------------------------------------------------------------------
# # Handle replicates
# #------------------------------------------------------------------------------
# #used to define the replicate structure
# _reps = {}
# for run in config['runs'].keys():
#     r = config['runs'][run]
#     tmp = []
#     for (rep, i) in enumerate(range(0, len(r), 2)):
#         if r[i]: tmp.append("rep%s" % str(rep+1))
#     _reps[run] = tmp
# #print(_reps)

#NOTE: Template class allows for _ in the variable names, we want to DISALLOW
#that for replicates
#ref: http://stackoverflow.com/questions/2326757/string-templates-in-python-what-are-legal-characters

class RepTemplate(Template):
    idpattern = r'[a-z][a-z0-9]*'

# #THIS helper fn is used in several of the modules peaks, ceas, frips
# #Instead of an expand, we need this fn to create the CORRECT input-list
# def _getRepInput(temp, suffix=""):
#     """generalized input fn to get the replicate files
#     CALLER passes in temp: a python string template that has the var runRep
#     e.g. analysis/ceas/$runRep/$runRep_DHS_stats.txt
#     Return: list of the string template filled with the correct runRep names
#     """
#     #print(temp)
#     s = RepTemplate(temp)
#     ls = []
#     for run in config['runs'].keys():
#         for rep in _reps[run]:
#             #GENERATE Run name: concat the run and rep name
#             runRep = "%s.%s" % (run, rep)
#             ls.append(s.substitute(runRep=runRep,))
#     #print(ls)
#     return ls


#################################################################
#######################======targets======#######################
#################################################################

def all_targets(wildcards):
    ls = []
    for sample in config['samples']:
        ls.append("analysis/%s/align/%s.sorted.bam" % (sample,sample))
        ls.append("analysis/%s/distribution/fragment_files/%s_fragment.txt" % (sample,sample))
        ls.append("analysis/%s/distribution/plots/%s_distribution.png" % (sample,sample))
        ls.append("analysis/%s/filtered/%s.sub50.sorted.bam" % (sample,sample))
        ls.append("analysis/%s/filtered/%s.sub100.sorted.bam" % (sample,sample))
        ls.append("analysis/%s/filtered/%s.sub150.sorted.bam" % (sample,sample))
        ls.append("analysis/%s/filtered/%s.sub200.sorted.bam" % (sample,sample))
        ls.append("analysis/%s/filtered/%s.raw.sorted.bam" % (sample,sample))
        for cutoff in ['raw','sub50','sub100','sub150','sub200']:
            ls.append("analysis/%s/peaks/%s.%s/%s.%s_peaks.narrowPeak" % (sample,sample,cutoff,sample,cutoff))
            ls.append("analysis/%s/peaks/%s.%s/%s.%s_peaks.bed" % (sample,sample,cutoff,sample,cutoff))
            ls.append("analysis/%s/peaks/%s.%s/%s.%s_summits.bed" % (sample,sample,cutoff,sample,cutoff))
            ls.append("analysis/%s/peaks/%s.%s/%s.%s_sorted_summits.bed" % (sample,sample,cutoff,sample,cutoff))
            ls.append("analysis/%s/peaks/%s.%s/%s.%s_sorted_5k_summits.bed" % (sample,sample,cutoff,sample,cutoff))
            ls.append("analysis/%s/peaks/%s.%s/%s.%s_peaks.xls" % (sample,sample,cutoff,sample,cutoff))
            ls.append("analysis/%s/peaks/%s.%s/%s.%s_treat_pileup.bdg" % (sample,sample,cutoff,sample,cutoff))
            ls.append("analysis/%s/peaks/%s.%s/%s.%s_control_lambda.bdg" % (sample,sample,cutoff,sample,cutoff))
            ls.append("analysis/%s/peaks/10FoldChange/%s.%s_peaks.narrowPeak" % (sample,sample,cutoff))
            ls.append("analysis/%s/peaks/20FoldChange/%s.%s_peaks.narrowPeak" % (sample,sample,cutoff))
            ls.append("analysis/%s/motif/%s.%s/" % (sample,sample,cutoff))
            ls.append("analysis/%s/motif/%s.%s/results" % (sample,sample,cutoff))
            ls.append("analysis/%s/motif/%s.%s/results/homerResults.html" % (sample,sample,cutoff))
            ls.append("analysis/%s/conservation/%s.%s/%s.%s_conserv.png" % (sample,sample,cutoff,sample,cutoff))
            ls.append("analysis/%s/DHS/%s.%s/%s.%s_DHS_stats.txt" % (sample,sample,cutoff,sample,cutoff))
    return ls


def getFastq(wildcards):
    return config["samples"][wildcards.sample]

def checkBAMPE(wildcards):
    first = config["samples"][wildcards.sample]
    ret = "-f BAMPE" if len(first) == 2 else ""
    return ret

def _createEmptyMotif(motif_html):
    """When the _sorted_5k_summits.bed has too few peaks, or is empty,
    we still want to create an emtpy homerResult.html
    INPUT: output paths of these files
    """
    #CHECK for dir existence:
    _path = "/".join(motif_html.split("/")[:-1])
    if not os.path.exists(_path):
        os.makedirs(_path, exist_ok=True)
    #Create an empty mdseqpos_index.html
    subprocess.call(['touch', motif_html])


#################################################################
#####################======parameters======######################
#################################################################

_threads = 8
_logfile = "analysis/logfile.log"
_macs_fdr="0.01"
_macs_keepdup="1"
_macs_extsize="146"
#_numPngs is used in conservation_plot rule to see how many pngs to expect
#note: the rule plots 3 runs per png, so for example, 12 runs results in 4 pngs
_nPerPlot = 3
_numPngs = math.ceil(len(config['samples'].keys())/float(_nPerPlot))
_nPngs = [n+1 for n in range(_numPngs)]


#################################################################
########################======start======########################
#################################################################


rule target_all:
    input:
        all_targets

#################################################################
########################======align======########################
#################################################################

rule align:
    input:
        getFastq
    output:
        temp("analysis/{sample}/align/{sample}.bam")
    params:
        index=config['bwa_index'],
    threads: _threads
    message: "ALIGN: Running BWA mem for alignment"
    log: _logfile
    shell:
        "bwa mem -t {threads} {params.index} {input} | samtools view -Sb - > {output} 2>>{log}"

rule sortBams:
    """General sort rule--take a bam {filename}.bam and 
    output {filename}.sorted.bam"""
    input:
        "analysis/{sample}/align/{sample}.bam"
        # getBam
    output:
        "analysis/{sample}/align/{sample}.sorted.bam",
        #"analysis/align/{sample}/{sample}.sorted.bam.bai"
    message: "ALIGN: sort bam file"
    log: _logfile
    threads: _threads
    shell:
        "sambamba sort {input} -o {output} -t {threads} 2>>{log}"

#################################################################
####################======distribution======#####################
#################################################################

rule get_fragments_length:
    input:
        "analysis/{sample}/align/{sample}.sorted.bam",
    output:
        "analysis/{sample}/distribution/fragment_files/{sample}_fragment.txt"
    message: "DISTRIBUTION: get frgament length"
    log: _logfile
    threads: _threads
    shell:
        "samtools view -q 1 -@ {threads} {input} | cut -f9 > {output} 2>>{log}"

rule plot_fragments_density:
    input:
        "analysis/{sample}/distribution/fragment_files/{sample}_fragment.txt"
    output:
        "analysis/{sample}/distribution/plots/{sample}_distribution.png"
    message: "DISTRIBUTION: get frgament length"
    log: _logfile
    shell:
        "Rscript ./FilterWorkFlow/scripts/density.r {input} {output}"

#################################################################
#######################======filter======########################
#################################################################

rule filter_bam_files_cutoff_50:
    input:
        "analysis/{sample}/align/{sample}.sorted.bam"
    output:
        "analysis/{sample}/filtered/{sample}.sub50.sorted.bam"
    message: "FILTERING: cutoff = 50"
    threads: _threads
    log: _logfile
    shell:
        "samtools view -h -@ {threads} {input} | awk '($9 <= 50 && $9 >= -50) || $1 ~ /^@/' | samtools view -bS -@ {threads} - > {output}"

rule filter_bam_files_cutoff_100:
    input:
        "analysis/{sample}/align/{sample}.sorted.bam"
    output:
        "analysis/{sample}/filtered/{sample}.sub100.sorted.bam"
    message: "FILTERING: cutoff = 100"
    threads: _threads
    log: _logfile
    shell:
        "samtools view -h -@ {threads} {input} | awk '($9 <= 100 && $9 >= -100) || $1 ~ /^@/' | samtools view -bS -@ {threads} - > {output}"

rule filter_bam_files_cutoff_150:
    input:
        "analysis/{sample}/align/{sample}.sorted.bam"
    output:
        "analysis/{sample}/filtered/{sample}.sub150.sorted.bam"
    message: "FILTERING: cutoff = 150"
    threads: _threads
    log: _logfile
    shell:
        "samtools view -h -@ {threads} {input} | awk '($9 <= 150 && $9 >= -150) || $1 ~ /^@/' | samtools view -bS -@ {threads} - > {output}"

rule filter_bam_files_cutoff_200:
    input:
        "analysis/{sample}/align/{sample}.sorted.bam"
    output:
        "analysis/{sample}/filtered/{sample}.sub200.sorted.bam"
    message: "FILTERING: cutoff = 200"
    threads: _threads
    log: _logfile
    shell:
        "samtools view -h -@ {threads} {input} | awk '($9 <= 200 && $9 >= -200) || $1 ~ /^@/' | samtools view -bS -@ {threads} - > {output}"

rule mv_raw_sorted_bam:
    input:
        "analysis/{sample}/align/{sample}.sorted.bam"
    output:
        "analysis/{sample}/filtered/{sample}.raw.sorted.bam"
    params:
        abspath=lambda wildcards, input: os.path.abspath(str(input))
    message: "soft link without filtering file to filtered folder"
    log: _logfile
    shell:
        "ln -s {params.abspath} {output}"

#################################################################
########################======peaks======########################
#################################################################

rule peaks_calling:
    input:
        "analysis/{sample}/filtered/{filename}.sorted.bam"
    output:
        "analysis/{sample}/peaks/{filename}/{filename}_peaks.narrowPeak",
        "analysis/{sample}/peaks/{filename}/{filename}_summits.bed",
        "analysis/{sample}/peaks/{filename}/{filename}_peaks.xls",
        "analysis/{sample}/peaks/{filename}/{filename}_treat_pileup.bdg",
        "analysis/{sample}/peaks/{filename}/{filename}_control_lambda.bdg",
    params:
        fdr=_macs_fdr,
        keepdup=_macs_keepdup,
        extsize=_macs_extsize,
        genome_size=config['genome_size'],
        outdir="analysis/{sample}/peaks/{filename}/",
        name="{filename}",
        #handle PE alignments--need to add -f BAMPE to macs2 callpeaks
        BAMPE = lambda wildcards: checkBAMPE(wildcards),
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        # treatment = lambda wildcards, input: [" -t %s" % i for i in input.treat] if input.treat else "",
        # control = lambda wildcards, input: [" -c %s" % i for i in input.cont] if input.cont else "",
    message: "PEAKS: calling peaks with macs2"
    log:_logfile
    shell:
       "{params.pypath} {config[macs2_path]} callpeak -B -q {params.fdr} --keep-dup {params.keepdup} -g {params.genome_size} {params.BAMPE} --extsize {params.extsize} --nomodel -t {input} --outdir {params.outdir} -n {params.name} 2>>{log}"

rule unsort_peaks_to_bed:
    input:
        "analysis/{sample}/peaks/{filename}/{filename}_peaks.narrowPeak"
    output:
        "analysis/{sample}/peaks/{filename}/{filename}_peaks.bed"
    message: "PEAKS: Converting unsorted peak file to bed file"
    log:_logfile
    shell:
        "cut -f1,2,3,4,9 {input} > {output} 2>>{log}"

rule sort_summits:
    input:
        "analysis/{sample}/peaks/{filename}/{filename}_summits.bed"
    output:
        "analysis/{sample}/peaks/{filename}/{filename}_sorted_summits.bed"
    message: "PEAKS: sorting the summits bed by score"
    log:_logfile
    shell:
        "sort -r -n -k 5 {input} > {output} 2>>{log}"

rule get_top_5k_summits:
    input:
        "analysis/{sample}/peaks/{filename}/{filename}_sorted_summits.bed"
    output:
        "analysis/{sample}/peaks/{filename}/{filename}_sorted_5k_summits.bed"
    message: "PEAKS: get top 5k summits bed by score"
    log:_logfile
    shell:
        "head -n 5000 {input} > {output} 2>>{log}"

rule fold_10_peaks:
    input:
        "analysis/{sample}/peaks/{filename}/{filename}_peaks.narrowPeak"
    output:
        "analysis/{sample}/peaks/10FoldChange/{filename}_peaks.narrowPeak"
    message: "PEAKS: Extract 10 fold change peaks"
    log:_logfile
    shell:
        "awk '$7 > 10' {input} > {output}"

rule fold_20_peaks:
    input:
        "analysis/{sample}/peaks/{filename}/{filename}_peaks.narrowPeak"
    output:
        "analysis/{sample}/peaks/20FoldChange/{filename}_peaks.narrowPeak"
    message: "PEAKS: Extract 20 fold change peaks"
    log:_logfile
    shell:
        "awk '$7 > 20' {input} > {output}"

#################################################################
########################======motif======########################
#################################################################

rule find_motif:
    input:
        bed="analysis/{sample}/peaks/{filename}/{filename}_sorted_5k_summits.bed"
    output:
        path="analysis/{sample}/motif/{filename}/",
        results="analysis/{sample}/motif/{filename}/results",
        html="analysis/{sample}/motif/{filename}/results/homerResults.html",
    params:
        genome=config['motif_path'],
        size=600,
    message: "MOTIF: calling HOMER on top 5k summits"
    threads:_threads
    log: _logfile
    run:
        #check to see if _sorted_5k_summits.bed is valid
        wc = str(subprocess.check_output(['wc', '-l', input.bed]))
        #this returns a byte string--we need to eliminate the b' ... '
        #and then convert the first elm to int
        wc = int(wc[2:-1].split()[0])
        if wc >= _minPeaks:
            #PASS- run motif scan
            shell("findMotifsGenome.pl {input} {params.genome} {output.results} -size {params.size} -p {threads} -mask -seqlogo -preparsedDir {output.results} >>{log} 2>&1")
        else:
            #FAIL - create empty outputs
            _createEmptyMotif(output.html)

#################################################################
####################======conservation======#####################
#################################################################

rule conservation:
    """generate conservation plots"""
    input:
        "analysis/{sample}/peaks/{filename}/{filename}_sorted_summits.bed"
    output:
        png="analysis/{sample}/conservation/{filename}/{filename}_conserv.png",
        thumb="analysis/{sample}/conservation/{filename}/{filename}_conserv_thumb.png",
        r="analysis/{sample}/conservation/{filename}/{filename}_conserv.R",
        score="analysis/{sample}/conservation/{filename}/{filename}_conserv.txt",
    params:
        db=config['conservation'],
        width=4000,
        #run = lambda wildcards: wildcards.run,
        # run="{filename}",
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
    message: "CONSERVATION: calling conservation script"
    log: _logfile
    shell:
        "{params.pypath} {config[python2]} ./FilterWorkFlow/scripts/conservation_plot.py -t Conservation_at_summits -d {params.db} -o analysis/{wildcards.sample}/conservation/{wildcards.filename}/{wildcards.filename}_conserv -l Peak_summits {input} -w {params.width} > {output.score} 2>>{log}"


#################################################################
#########################======DHS======#########################
#################################################################

rule DHS_intersectDHS:
    """Intersect PEAKS with DHS regions"""
    input:
        'analysis/{sample}/peaks/{filename}/{filename}_peaks.bed'
    output:
        'analysis/{sample}/DHS/{filename}/{filename}_DHS_peaks.bed'
    params:
        #check for config['DHS'] defined, otherwise, use null
        dhs=config['DHS'] if config['DHS'] else "/dev/null"
    message: "DHS: intersect PEAKS with DHS regions"
    log: _logfile
    shell:
        "intersectBed -wa -u -a {input} -b {params.dhs} > {output} 2>>{log}"

rule DHS_stat:
    """collect DHS stats"""
    input:
        n='analysis/{sample}/peaks/{filename}/{filename}_peaks.bed',
        dhs='analysis/{sample}/DHS/{filename}/{filename}_DHS_peaks.bed'
    output:
        'analysis/{sample}/DHS/{filename}/{filename}_DHS_stats.txt'
    message: "DHS: collecting stats"
    log: _logfile
    shell:
        "wc -l {input.n} {input.dhs} > {output} 2>>{log}"














