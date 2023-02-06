import pysam
import os
import numpy as np

def write_fastq (read, fastqfile_fo):
    fastqfile_fo.write("@%s\n" % read.qname)
    fastqfile_fo.write("%s\n" % read.seq)
    fastqfile_fo.write("+\n")
    fastqfile_fo.write("%s\n" % read.qual)
    fastqfile_fo.flush()

def write_fastq_original (read, fastqfile_fo):
    if not read.is_reverse:
        write_fastq(read, fastqfile_fo)
    else:
        fastqfile_fo.write("@%s\n" % read.qname)
        fastqfile_fo.write("%s\n" % read.get_forward_sequence())
        fastqfile_fo.write("+\n")
        for_qualities = []
        rev_qualities = read.qual
        for i in range(0, len(rev_qualities)):
            for_qualities.append(rev_qualities[len(rev_qualities) - i - 1])
        fastqfile_fo.write("%s\n" % "".join(for_qualities))
        fastqfile_fo.flush()

def assemble_reads_lamassemble (path_to_fastq, path_to_lamassemble, path_to_lamassemble_train_file, assembly_name, savepath_for_assembly_fasta):
    os.system(path_to_lamassemble + " -P 8 --all -n " + assembly_name + " " + path_to_lamassemble_train_file + " " + path_to_fastq + " > " + savepath_for_assembly_fasta)

def assemble_reads_lamassemble_s3 (path_to_fastq, path_to_lamassemble, path_to_lamassemble_train_file, assembly_name, savepath_for_assembly_fasta):
    os.system(path_to_lamassemble + " -P 16 -s 3 -n " + assembly_name + " " + path_to_lamassemble_train_file + " " + path_to_fastq + " > " + savepath_for_assembly_fasta)

def assemble_reads_lamassemble_s (path_to_fastq, path_to_lamassemble, path_to_lamassemble_train_file, assembly_name, savepath_for_assembly_fasta, s):
    if s == "1":
        os.system(path_to_lamassemble + " -P 16 --all -n " + assembly_name + " " + path_to_lamassemble_train_file + " " + path_to_fastq + " > " + savepath_for_assembly_fasta)
    else:
        os.system(path_to_lamassemble + " -P 16 -s " + s + " -n " + assembly_name + " " + path_to_lamassemble_train_file + " " + path_to_fastq + " > " + savepath_for_assembly_fasta)

def collect_reads (read_names, chrom, pos, search_area, alignment_file):

    reads = {}

    for read in alignment_file.fetch(chrom, np.amax([0, pos - search_area]), pos + search_area):
        if read.query_name in read_names and not (read.query_name in reads):
            reads[read.query_name] = read

    return reads

def read_single_fasta (fasta_file_path):
    fasta_name = ""
    fasta_sequence = ""
    for line in open(fasta_file_path):
        if line.startswith(">"):
            if fasta_name == "":
                fasta_name = line[1:].rstrip()
            else:
                break
        else:
            fasta_sequence += line.rstrip()
    return fasta_name, fasta_sequence
