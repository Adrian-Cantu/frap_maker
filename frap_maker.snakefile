import os
import sys
import socket

hostname = socket.gethostname()

sys.stderr.write(f"Running on: {hostname}\n")

#configfile: 'frap_maker.json'


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
        expand(os.path.join(SMALT_OUT, "{sample}_good_out_vs_"+ DATABASE +".sam"), sample=SAMPLES),
        expand(os.path.join(SMALT_OUT, "vs_"+ DATABASE +"_hits.{sample}_good_out.tab"),sample=SAMPLES),
#        os.path.join(SMALT_OUT,'all_normalized.txt'),
#        os.path.join(SMALT_OUT,'all_hits.txt'),
        os.path.join(SMALT_OUT,'all_normalized_per_million.txt')

rule smalt:
    input:
        fasta=os.path.join(FASTA, "{sample}_good_out." + FQEXTN )
    output:
        os.path.join(SMALT_OUT, "{sample}_good_out_vs_" + DATABASE + ".sam")
    params:
        outfile=os.path.join(SMALT_OUT, "{sample}_good_out_vs_"+ DATABASE +".sam"),
        database=os.path.join(DB_DIR, DATABASE)
    threads: THREADS
    shell:
        "{config[executables][smalt]}  map -n {threads} -f sam -y 1 -o {params.outfile} {params.database} {input.fasta} "

rule gen_tsv:
    input:
        sam=os.path.join(SMALT_OUT, "{sample}_good_out_vs_"+ DATABASE +".sam")
    output:
        os.path.join(SMALT_OUT, "vs_"+ DATABASE +"_hits.{sample}_good_out.tab")
    params:
        outfile=os.path.join(SMALT_OUT, "vs_"+ DATABASE +"_hits.{sample}_good_out.tab")
    threads: 1
    shell:
        'grep -v ^@ {input.sam} | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c | sort -nr  | sed -e "s/^ *//" | tr " " "\\t"  > {params.outfile}'

rule do_tj:
    output:
        os.path.join(SMALT_OUT,'tj.txt')
    params:
        datadir = FASTA,
        outdir  = SMALT_OUT
    threads: 1
    shell:
        'perl do_tj.pl {params.datadir} {params.outdir}'

rule smalt_norm:
    input:
        tj = os.path.join(SMALT_OUT,'tj.txt'),
        files = expand(os.path.join(SMALT_OUT, "vs_"+ DATABASE +"_hits.{sample}_good_out.tab"),sample=SAMPLES)
    output:
#        f1 = os.path.join(SMALT_OUT,'all_normalized.txt'),
#        f2 = os.path.join(SMALT_OUT,'all_hits.txt'),
        f3 = os.path.join(SMALT_OUT,'all_normalized_per_million.txt')
    params:
        database=os.path.join(DB_DIR, DATABASE)
    threads: 1
    shell:
        "perl frap_normalization_better.pl -t {input.tj} -m -l 50000 -f {params.database} {input.files} > {output.f3}"
        


#rule frap_normal:
#    input:
        
