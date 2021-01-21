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
        expand(os.path.join(SMALT_OUT, "{sample}_good_out_vs_"+ DATABASE +".sam"), sample=SAMPLES)


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

