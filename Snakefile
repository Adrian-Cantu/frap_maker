import os
import sys
import socket

hostname = socket.gethostname()

sys.stderr.write(f"Running on: {hostname}\n")

configfile: 'sample.json'


FASTA       = config['directories']['data']
SMALT_OUT   = config['directories']['smalt_out']
THREADS     = config['threads']
DB_DIR      = config['directories']['db_dir']
DATABASE    = config['directories']['db']

if hostname.startswith('node'):
    hostname = 'anth'
else:
    hostname = 'tatabox'

if f"{hostname}_executables" not in config:
    sys.stderr.write(f"No executables defined for {hostname}_executables in process_metagenomes.json\n")
    sys.exit(-1)

config['executables'] = config[f"{hostname}_executables"]

# A Snakemake regular expression matching the forward mate FASTQ files.
# the comma after SAMPLES is important!
SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(FASTA, '{sample}_good_out.{extn}'))
if len(SAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
    sys.exit(0)
if len(set(EXTENSIONS)) != 1:
    sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
    sys.stderr.write("\n\t".join(set(EXTENSIONS)))
    sys.stderr.write("\nWe don't know how to handle these\n")
    sys.exit(0)

FQEXTN = EXTENSIONS[0]
PATTERN_R1 = '{sample}_good_out.' + FQEXTN

print(f"Samples are {SAMPLES}")

# this is the samples we have processed!

rule all:
    input:
        expand(os.path.join(SMALT_OUT, f"{{sample}}_good_out_vs_{DATABASE}.sam"), sample=SAMPLES),
        expand(os.path.join(SMALT_OUT, f"vs_{DATABASE}_hits.{{sample}}_good_out.tab"),sample=SAMPLES),
        os.path.join(SMALT_OUT,'all_normalized.txt'),
        os.path.join(SMALT_OUT,'all_hits.txt'),
        os.path.join(SMALT_OUT,'all_normalized_per_million.txt')

rule smalt_map:
    input:
        fasta=os.path.join(FASTA, f"{{sample}}_good_out.{FQEXTN}" )
    output:
        outfile=os.path.join(SMALT_OUT, f"{{sample}}_good_out_vs_{DATABASE}.sam")
    params:
        database=os.path.join(DB_DIR, DATABASE),
        smalt=config['executables']['smalt'],
        logs="logs/smalt_map.{sample}.log"
    conda: 
        "environment.yml"
    threads: THREADS
    log: "logs/smalt_map.{sample}.log"
    shell:
        "{params.smalt} map -n {threads} -f sam -y 1 -o {output.outfile} {params.database} {input.fasta} 2> {params.logs} "

rule gen_tsv:
    input:
        sam=os.path.join(SMALT_OUT, f"{{sample}}_good_out_vs_{DATABASE}.sam")
    output:
        outfile=os.path.join(SMALT_OUT, f"vs_{DATABASE}_hits.{{sample}}_good_out.tab")
    conda:
        "environment.yml"
    threads: 1
    log: "logs/gen_tsv.{sample}.log"
    shell:
        'grep -v ^@ {input.sam} | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c | sort -nr  | sed -e "s/^ *//" | tr " " "\\t"  > {output.outfile} 2> {log}'

rule do_tj:
    output:
        os.path.join(SMALT_OUT,'tj.txt')
    params:
        datadir = FASTA,
        outdir= lambda wildcards, output: os.path.splitext(output[0])[0]
#        outdir  = SMALT_OUT
    conda:
        "environment.yml"
    threads: 1
    log: "logs/do_tj.log"
    shell:
        'perl do_tj.pl {params.datadir} {params.outdir} 2> {log}'

rule smalt_norm:
    input:
        tj = os.path.join(SMALT_OUT,'tj.txt'),
        files = expand(os.path.join(SMALT_OUT, f"vs_{DATABASE}_hits.{{sample}}_good_out.tab"),sample=SAMPLES)
    output:
        f1 = os.path.join(SMALT_OUT,'all_normalized.txt'),
        f2 = os.path.join(SMALT_OUT,'all_hits.txt'),
        f3 = os.path.join(SMALT_OUT,'all_normalized_per_million.txt')
    params:
        database=os.path.join(DB_DIR,f"{DATABASE}.fasta" )
    conda:
        "environment.yml"
    log:
        "logs/smalt_norm.log"
    threads: 1
    shell:
        """
        perl frap_normalization_better.pl -t {input.tj} -n -l 50000 -f {params.database} {input.files} > {output.f1} 2> {log}
        perl frap_normalization_better.pl -t {input.tj} -h -f {params.database} {input.files} > {output.f2} 2> {log}
        perl frap_normalization_better.pl -t {input.tj} -m -l 50000 -f {params.database} {input.files} > {output.f3} 2> {log}
        """
        


#rule frap_normal:
#    input:
        
