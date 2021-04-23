import os
import sys
import socket

hostname = socket.gethostname()

sys.stderr.write(f"Running on: {hostname}\n")

configfile: 'sample.json'


FASTA       = config['directories']['data']
SMALT_OUT   = config['directories']['smalt_out']
THREADS     = config['threads']
IDENT       = config['identity']
DB_DIR      = config['directories']['db_dir']
DATABASE = config['directories']['db']
DATABASE_FILE_BASE   = os.path.splitext(DATABASE)[0]

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

#this part defines how fastas are scanned
#SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(FASTA, '{sample}_good_out.{extn}'))
SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(FASTA, '{sample}.sra_1.{extn}'))

if len(SAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {FASTA}. Is this the right read dir?\n")
    sys.exit(0)
if len(set(EXTENSIONS)) != 1:
    sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
    sys.stderr.write("\n\t".join(set(EXTENSIONS)))
    sys.stderr.write("\nWe don't know how to handle these\n")
    sys.exit(0)

FQEXTN = EXTENSIONS[0]
#PATTERN_R1 = '{sample}_good_out.' + FQEXTN

print(f"Samples are {SAMPLES}")

# this is the samples we have processed!

rule all:
    input:
        expand(os.path.join(SMALT_OUT, f"{{sample}}_good_out_vs_{DATABASE_FILE_BASE}.sam"), sample=SAMPLES),
        expand(os.path.join(SMALT_OUT, f"vs_{DATABASE_FILE_BASE}_hits.{{sample}}_good_out.tab"),sample=SAMPLES),
        os.path.join(SMALT_OUT,'all_normalized.txt'),
        os.path.join(SMALT_OUT,'all_hits.txt'),
        os.path.join(SMALT_OUT,'all_normalized_per_million.txt')
#        expand(os.path.join(FASTA, f"{{sample}}_good_out.fasta" ), sample=SAMPLES)

rule prinseq:
    input:
        infasta = os.path.join(FASTA, f"{{sample}}.sra_1.{FQEXTN}" )
    output:
        outfile=os.path.join(FASTA, f"{{sample}}_good_out.fasta" )
    params:
        #sample=os.path.join(FASTA,"{sample}")
        sample= lambda wildcards, output: os.path.join(os.path.dirname(output[0]),f"{wildcards.sample}")
    conda:
        "environment.yml"
    threads: THREADS
    log: "logs/prinseq.{sample}.log"
    shell:
        "prinseq++ -fastq {input.infasta} -min_len 50 -min_qual_mean 30 -lc_entropy=0.5 -out_format 1 -out_name {params.sample} -threads {threads} 2> {log}"

rule smalt_index:
    input:
        os.path.join(DB_DIR, DATABASE)
    output:
        sma=os.path.join(DB_DIR,f"{DATABASE_FILE_BASE}.sma"),
        smi=os.path.join(DB_DIR,f"{DATABASE_FILE_BASE}.smi")
    params:
        base=lambda wildcards, output: os.path.splitext(output[0])[0]
    conda:
        "environment.yml"
    threads: 1
    log: "logs/smalt_index.log"
    shell:
        "smalt index -k 10 -s 5 {params.base} {input} 2> {log}"


rule smalt_map:
    input:
        fasta=os.path.join(FASTA, f"{{sample}}_good_out.fasta" ),
        sma=os.path.join(DB_DIR,f"{DATABASE_FILE_BASE}.sma"),
        smi=os.path.join(DB_DIR,f"{DATABASE_FILE_BASE}.smi")
    output:
        outfile=os.path.join(SMALT_OUT, f"{{sample}}_good_out_vs_{DATABASE_FILE_BASE}.sam")
    params:
        database=lambda wildcards, input: os.path.splitext(f"{input.sma}")[0],
        smalt=config['executables']['smalt'],
        logs="logs/smalt_map.{sample}.log",
        idd=IDENT
    conda: 
        "environment.yml"
    threads: THREADS
    log: "logs/smalt_map.{sample}.log"
    shell:
        "smalt map -n {threads} -f sam -y {params.idd} -o {output.outfile} {params.database} {input.fasta} 2> {params.logs} "

rule gen_tsv:
    input:
        sam=os.path.join(SMALT_OUT, f"{{sample}}_good_out_vs_{DATABASE_FILE_BASE}.sam")
    output:
        outfile=os.path.join(SMALT_OUT, f"vs_{DATABASE_FILE_BASE}_hits.{{sample}}_good_out.tab")
    conda:
        "environment.yml"
    threads: 1
    log: "logs/gen_tsv.{sample}.log"
    shell:
        'grep -v ^@ {input.sam} | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c | sort -nr  | sed -e "s/^ *//" | tr " " "\\t"  > {output.outfile} 2> {log} && [[ -s {output.outfile} ]]'

rule do_tj:
    input:
        expand(os.path.join(FASTA, f"{{sample}}_good_out.fasta" ), sample=SAMPLES)
    output:
        os.path.join(SMALT_OUT,'tj.txt')
    params:
        datadir = lambda wildcards, input: os.path.dirname(input[0]),
        outdir= lambda wildcards, output: os.path.dirname(output[0])
#        outdir  = SMALT_OUT
    conda:
        "environment.yml"
    threads: 1
    log: "logs/do_tj.log"
    shell:
        'perl do_tj.pl {params.datadir} {params.outdir} 2> {log} && [[ -s {output} ]]'

rule smalt_norm:
    input:
        tj = os.path.join(SMALT_OUT,'tj.txt'),
        files = expand(os.path.join(SMALT_OUT, f"vs_{DATABASE_FILE_BASE}_hits.{{sample}}_good_out.tab"),sample=SAMPLES)
    output:
        f1 = os.path.join(SMALT_OUT,'all_normalized.txt'),
        f2 = os.path.join(SMALT_OUT,'all_hits.txt'),
        f3 = os.path.join(SMALT_OUT,'all_normalized_per_million.txt')
    params:
        database=os.path.join(DB_DIR,f"{DATABASE}" )
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
        
