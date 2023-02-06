import numpy as np
import re
import vcf
import read_methods
import mappy as mp
import os
from time import time_ns
import pysam

# gather reads that span an insertion (either as one mapping or as primary and supplementary mappings)
# samfile : pysam samfile of the sample
# sv_record : vcf record of this sv
# sv_breakpoint : string either start or end depending which breakpoint of the sv should the analyzed
# debug : bool if true write debug messages and bed files
def gather_germline_insertion_evidence (samfile, sv_record, sv_breakpoint, debug):

    time_total_start = time_ns()
    time_single_spanning = 0
    time_SA_spanning = 0
    A = 0
    B = 0

    sv_chrom = sv_record.CHROM
    sv_pos = sv_record.POS
    if sv_breakpoint == "END":
        if sv_record.INFO["SVTYPE"] == "BND":
            sv_chrom = sv_record.ALT[0].chr
            sv_pos = sv_record.ALT[0].pos
        else:
            sv_chrom = sv_record.INFO["CHR2"]
            sv_pos = sv_record.INFO["END"]
    search_area = 3000

    if debug:
        bed_out = open("tmp.bed", "w")

    reads_of_interest = {}

    if debug:
        print("Looking at " + sv_chrom + ":" + str(np.amax([0, sv_pos - search_area])) + "-" + str(sv_pos + search_area))

    for read in samfile.fetch(sv_chrom, np.amax([0, sv_pos - search_area]), sv_pos + search_area):

        # Reads spanning the breakpoint as one primary mapping, at least 200 bp left and right of the breakpoint
        if sv_pos - read.reference_start > 200 and read.reference_end - sv_pos > 200:

            has_long_insertion = False

            cigar_pos = read.reference_start

            # crawl through the read using cigar operations
            for cigar in read.cigartuples:
                # add M and D to current position
                if cigar[0] in [0, 2]:
                    cigar_pos += cigar[1]
                # if insertion longer than 200 bp
                if cigar[0] == 1 and cigar[1] > 200:
                    # if distance between insertion position and breakpoint is shorter than 200 bp
                    if np.abs(cigar_pos - sv_pos) < 200:
                        has_long_insertion = True

            if has_long_insertion:
                if read.query_name in reads_of_interest:
                    reads_of_interest[read.query_name].append(read)
                else:
                    reads_of_interest[read.query_name] = [read]

        # Reads spanning the SV as primary and supplementary mapping
        if read.has_tag("SA"):

            read_orientation = "-"
            if read.is_forward:
                read_orientation = "+"

            # read does not span the breakpoint
            if not (read.reference_start + 200 < sv_pos and read.reference_end > sv_pos + 200):

                # check for each read + SA combination if the combined interval spans the breakpoint and has the same orientation
                for SA in read.get_tag("SA").split(";"):

                    if SA != "":

                        # get SA informations
                        SA_chrom = SA.split(",")[0]
                        SA_start = int(SA.split(",")[1])
                        SA_cigar = SA.split(",")[3]
                        SA_end = SA_start
                        SA_orientation = SA.split(",")[2]
                        if read.query_name == "ac9ca1f5-288f-4f8c-9988-68f3dfca618f":
                            print(SA)

                        # infer alignment end position by crawling through the cigar string
                        for op in re.findall("[0-9]+[A-Z]{1}", SA_cigar):
                            if op[-1] in ["M", "D"]:
                                SA_end += int(op[:-1])

                        # check if SA is on the same chromosome
                        if SA_chrom == sv_chrom and SA_orientation == read_orientation:

                            a_bed_fo = open("tmp.a.bed", "w")
                            a_bed_fo.write("\t".join([read.reference_name, str(read.reference_start), str(read.reference_end)]) + "\n")
                            a_bed_fo.close()

                            b_bed_fo = open("tmp.b.bed", "w")
                            b_bed_fo.write("\t".join([SA_chrom, str(SA_start), str(SA_end)]) + "\n")
                            b_bed_fo.close()

                            os.system("bedtools closest -d -a tmp.a.bed -b tmp.b.bed > tmp.closest.bed")
                            distance = int(open("tmp.closest.bed", "r").read().split("\t")[-1])

                            if distance < 200:

                                os.system("cat tmp.a.bed tmp.b.bed | bedtools sort -i - | bedtools merge -d " + str(1 + distance) + " -i - > tmp.merge.bed")
                                merge_start = int(open("tmp.merge.bed", "r").read().split("\t")[1])
                                merge_end = int(open("tmp.merge.bed", "r").read().split("\t")[2])

                                # check if alignment spans the breakpoint
                                if merge_start + 200 < sv_pos and merge_end > sv_pos + 200:
                                    if read.query_name in reads_of_interest:
                                        reads_of_interest[read.query_name].append(read)
                                    else:
                                        reads_of_interest[read.query_name] = [read]

    if debug:
        for key in reads_of_interest:
            if len(reads_of_interest[key]) == 1:
                read = reads_of_interest[key][0]
                bed_out.write("\t".join(
                    [
                        read.reference_name,
                        str(read.reference_start),
                        str(read.reference_end),
                        key,
                        "1",
                        ".",
                        str(read.reference_start),
                        str(read.reference_end),
                        "0,0,255"
                    ]) + "\n")
            else:
                start = 1000000000
                end = 0
                for read in reads_of_interest[key]:
                    if read.reference_start < start:
                        start = read.reference_start
                    if read.reference_end > end:
                        end = read.reference_end
                bed_out.write("\t".join(
                    [
                        read.reference_name,
                        str(start),
                        str(end),
                        key,
                        "1",
                        ".",
                        str(start),
                        str(end),
                        "255,0,0"
                    ]) + "\n")

    if len(reads_of_interest) > 0:
        insertion_supporting_reads = []
        for key in reads_of_interest:
            insertion_supporting_reads.append(reads_of_interest[key][0])
        return insertion_supporting_reads

    else:
        return -1

# check if a sequences mapped using the aligner spans the position (chrom:pos) with at least padding on each site
def sequence_spans_position (chrom, pos, aligner, sequence, padding):
    sequence_spans_postion = False
    alignments_list = list(aligner.map(sequence))
    for aln in alignments_list:
        if aln.mapq > 0:
            if aln.ctg == chrom and aln.r_st + padding < pos and aln.r_en - padding > pos:
                sequence_spans_postion = True
    return sequence_spans_postion

def germline_filter (reference_fasta_file_path, prefiltered_vcf_file_path, bam_file_path, data_folder):

    ###
    # building folders
    # sv_fastq, sv_assemblies, sv_mappings
    # germline_fastq, germline_assemblies, germline_mappings
    ###

    if not os.path.exists(data_folder + "sv_fastq/"):
        os.system("mkdir " + data_folder + "sv_fastq/")
    if not os.path.exists(data_folder + "sv_assemblies/"):
        os.system("mkdir " + data_folder + "sv_assemblies/")
    if not os.path.exists(data_folder + "sv_mappings/"):
        os.system("mkdir " + data_folder + "sv_mappings/")
    if not os.path.exists(data_folder + "germline_fastq/"):
        os.system("mkdir " + data_folder + "germline_fastq/")
    if not os.path.exists(data_folder + "germline_assemblies/"):
        os.system("mkdir " + data_folder + "germline_assemblies/")
    if not os.path.exists(data_folder + "germline_mappings/"):
        os.system("mkdir " + data_folder + "germline_mappings/")

    ###
    # germline filter
    ###

    time_aligner_load_start_ns = time_ns()

    aligner = mp.Aligner(reference_fasta_file_path, preset="map-ont")
    if not aligner: raise Exception("ERROR: failed to load/build index")

    time_aligner_load_end_ns = time_ns()

    print("Loading aligner in " + str(round((time_aligner_load_end_ns - time_aligner_load_start_ns) / 1000000000, 2)) + " s")

    vcf_reader = vcf.Reader(open(data_folder + prefiltered_vcf_file_path, "r"))
    vcf_writer = vcf.Writer(open(data_folder + prefiltered_vcf_file_path.replace(".prefilter.vcf", ".filter.vcf"), "w"), template=vcf_reader)

    samfile = pysam.AlignmentFile(bam_file_path, "rb")

    for record in vcf_reader:
        if record.FILTER == []:
            #print("Starting analysis for SV record " + record.ID)

            ###
            # assemble sv consensus from reads using lamassemble
            ###

            time_sv_assembly_start_ns = time_ns()

            # collect sv reads in fastq
            sv_reads = read_methods.collect_reads(record.INFO["RNAMES"], record.CHROM, record.POS, 3000, samfile)
            fastq_fo = open(data_folder + "sv_fastq/" + record.ID + ".sv.fastq", "w")
            for key in sv_reads:
                read_methods.write_fastq_original(sv_reads[key], fastq_fo)
            fastq_fo.close()

            # assemble consensus
            read_methods.assemble_reads_lamassemble_s \
                    (
                    data_folder + "sv_fastq/" + record.ID + ".sv.fastq",
                    "lamassemble",
                    "/home/ubuntu/seq/Miltenyi_ONT/lamassemble/lamassemble/train/promethion.mat",
                    record.ID + ".sv",
                    data_folder + "sv_assemblies/" + record.ID + ".sv.fasta",
                    "2"
                )

            ###
            # check if sv assembly spans one breakpoint and indicates insertion
            ###

            assembly_sequence = list(mp.fastx_read(data_folder + "sv_assemblies/" + record.ID + ".sv.fasta", read_comment=False))[0][1]

            if record.INFO["SVTYPE"] == "BND":
                if np.any([
                    sequence_spans_position(record.CHROM, record.POS, aligner, assembly_sequence, 200),
                    sequence_spans_position(record.ALT[0].chr, record.ALT[0].pos, aligner, assembly_sequence, 200)
                ]):
                    record.FILTER.append("ASSEMBLY_SPAN")
            else:
                if np.any([
                    sequence_spans_position(record.CHROM, record.POS, aligner, assembly_sequence, 200),
                    sequence_spans_position(record.INFO["CHR2"], record.INFO["END"], aligner, assembly_sequence, 200)
                ]):
                    record.FILTER.append("ASSEMBLY_SPAN")

            time_sv_assembly_end_ns = time_ns()

            #print("\tBuilding SV assembly and checking for breakpoint span in " + str(round((time_sv_assembly_end_ns - time_sv_assembly_start_ns) / 1000000000, 2)) + " s")

            time_germline_analysis_start_ns = time_ns()

            ###
            # germline test
            ###
            #print("\tStart looking for evidence of germline insertion")
            for breakpoint in ["START", "END"]:
                reads = gather_germline_insertion_evidence(samfile, record, breakpoint, False)
                if reads != -1 and len(reads) <= 50:
                    #print("\tFound " + str(len(reads)) + " reads at " + breakpoint + " position")
                    fastq_fo = open(data_folder + "germline_fastq/" + record.ID + "." + breakpoint + ".germline_ins.fastq", "w")
                    for read in reads:
                        read_methods.write_fastq_original(read, fastq_fo)
                    fastq_fo.close()

                    # assemble consensus
                    overlap_cutoff = "2"
                    if len(reads) < 3:
                        overlap_cutoff = "1"
                    read_methods.assemble_reads_lamassemble_s \
                            (
                            data_folder + "germline_fastq/" + record.ID + "." + breakpoint + ".germline_ins.fastq",
                            "lamassemble",
                            "/home/ubuntu/seq/Miltenyi_ONT/lamassemble/lamassemble/train/promethion.mat",
                            record.ID + "." + breakpoint + ".germline_ins",
                            data_folder + "germline_assemblies/" + record.ID + "." + breakpoint + ".germline_ins.fasta",
                            overlap_cutoff
                        )

                    germline_assembly_sequence = list(mp.fastx_read(data_folder + "germline_assemblies/" + record.ID + "." + breakpoint + ".germline_ins.fasta",read_comment=False))[0][1]

                    if record.INFO["SVTYPE"] == "BND":
                        if np.any([
                            sequence_spans_position(record.CHROM, record.POS, aligner, germline_assembly_sequence, 200),
                            sequence_spans_position(record.ALT[0].chr, record.ALT[0].pos, aligner, germline_assembly_sequence, 200)
                        ]):
                            record.FILTER.append("GERMLINE_SPAN_" + breakpoint)
                    else:
                        if np.any([
                            sequence_spans_position(record.CHROM, record.POS, aligner, germline_assembly_sequence, 200),
                            sequence_spans_position(record.INFO["CHR2"], record.INFO["END"], aligner, germline_assembly_sequence, 200)
                        ]):
                            record.FILTER.append("GERMLINE_SPAN_" + breakpoint)
                if reads != -1 and len(reads) > 50:
                    record.FILTER.append("HIGH_COVERAGE")
            time_germline_analysis_end_ns = time_ns()

            #print("\tSpend " + str(round((time_germline_analysis_end_ns - time_germline_analysis_start_ns) / 1000000000, 2)) + " s gathering evidence for germline insertion and building germline assembly")

        vcf_writer.write_record(record)
