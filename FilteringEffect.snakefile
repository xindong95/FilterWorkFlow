#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import yaml

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
config = getRuns(config)
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

#THIS helper fn is used in several of the modules peaks, ceas, frips
#Instead of an expand, we need this fn to create the CORRECT input-list
def _getRepInput(temp, suffix=""):
    """generalized input fn to get the replicate files
    CALLER passes in temp: a python string template that has the var runRep
    e.g. analysis/ceas/$runRep/$runRep_DHS_stats.txt
    Return: list of the string template filled with the correct runRep names
    """
    #print(temp)
    s = RepTemplate(temp)
    ls = []
    for run in config['runs'].keys():
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            ls.append(s.substitute(runRep=runRep,))
    #print(ls)
    return ls

#------------------------------------------------------------------------------
# TARGETS
#------------------------------------------------------------------------------

def all_targets(wildcards):
    # _qdnaseq = config["cnv_analysis"]
    ls = []
    # #IMPORT all of the module targets
    # ls.extend(align_targets(wildcards))
    # ls.extend(peaks_targets(wildcards))
    # ls.extend(fastqc_targets(wildcards))
    # ls.extend(conservation_targets(wildcards))
    # ls.extend(ceas_targets(wildcards))
    # ls.extend(frips_targets(wildcards))
    # ls.extend(contamination_targets(wildcards))
    # ls.extend(mapmaker_targets(wildcards))
    # #ls.extend(bam_snapshots_targets(wildcards))
    # ls.extend(regulatory_targets(wildcards))

    # #Check to see if motif is enabled
    # if 'motif' in config:
    #     ls.extend(motif_targets(wildcards))

    # #HANDLE CNV/qdnaseq analysis
    # if _qdnaseq:
    #     #ls.extend(qdnaseq_targets(wildcards))

    #     #check for some inputs
    #     hasInput = False

    #     #HACK: for some reason, using the following line causes errors
    #     #for (run, ls) in config['runs'].items():
    #     #SO we call getRuns (from above) using a simplified config
    #     tmp_config = {'metasheet': config['metasheet']}
    #     runs = getRuns(tmp_config)['runs'].copy()
    #     for (run) in runs.keys():
    #         #NOTE: if i do this, this is an error!
    #         #ls = runs[run]
    #         if runs[run][1] or runs[run][3]:
    #             #these are the control sample indices
    #             hasInput = True
    #             break

    #     if hasInput:
    #         ls.extend(qdnaseq_targets(wildcards))

    # if "epicypher_analysis" in config and config["epicypher_analysis"]:
    #     ls.extend(epicypher_targets(wildcards))
        
    # ls.extend(report_targets(wildcards))
    return ls

# include:
def getFastq(wildcards):
    return config["samples"][wildcards.sample]

_bwa_threads = 8


rule target_all:
    all_targets

rule align:
    input:
        getFastq
    output:
        temp("analysis/align/{sample}/{sample}.bam")
    params:
        index=config['bwa_index'],
    threads: _bwa_threads
    message: "ALIGN: Running BWA mem for alignment"
    log: _logfile
    shell:
        "bwa mem -t {threads} {params.index} {input} | samtools view -Sb - > {output} 2>>{log}"

















