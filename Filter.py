import vcf
import numpy as np
import Misc
import pysam
import Read_methods
import os

def reformat_vcf (vcf_file_path, reformat_vcf_file_path):

    reformat_fo = open(reformat_vcf_file_path, "w")

    for line in open(vcf_file_path):
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                reformat_fo.write("##INFO=<ID=MQ_START,Number=1,Type=Float,Description=\"Mean mapping quality of supporting reads at SV start\">\n")
                reformat_fo.write("##INFO=<ID=MQ_END,Number=1,Type=Float,Description=\"Mean mapping quality of supporting reads at SV end\">\n")
                reformat_fo.write("##INFO=<ID=RMDUP_RNAMES,Number=.,Type=String,Description=\"Names of supporting reads without duplicates\">\n")
                reformat_fo.write("##INFO=<ID=RMDV,Number=1,Type=Integer,Description=\"Number of variant reads without duplicates\">\n")
                reformat_fo.write("##INFO=<ID=MISS_START,Number=1,Type=Float,Description=\"Number of mismatches divided by aligned sequence length for the reads at SV start\">\n")
                reformat_fo.write("##INFO=<ID=MISS_END,Number=1,Type=Float,Description=\"Number of mismatches divided by aligned sequence length for the reads at SV end\">\n")
                reformat_fo.write("##INFO=<ID=MISS_START_NORM,Number=1,Type=Float,Description=\"Number of mismatches divided by aligned sequence length for the non SV reads at SV start\">\n")
                reformat_fo.write("##INFO=<ID=MISS_END_NORM,Number=1,Type=Float,Description=\"Number of mismatches divided by aligned sequence length for the non SV reads at SV end\">\n")
                reformat_fo.write("##INFO=<ID=STD_START,Number=1,Type=Float,Description=\"Standard deviation of read end towards the SV breakpoint at the SV start\">\n")
                reformat_fo.write("##INFO=<ID=STD_END,Number=1,Type=Float,Description=\"Standard deviation of read end towards the SV breakpoint at the SV end\">\n")
                reformat_fo.write("##INFO=<ID=DIST_START,Number=1,Type=Integer,Description=\"Read start distance, only written for SVs with DV and RMDV of 2 from the same library\">\n")
                reformat_fo.write("##INFO=<ID=DIST_END,Number=1,Type=Integer,Description=\"Read end distance, only written for SVs with DV and RMDV of 2 from the same library\">\n")
                reformat_fo.write("##INFO=<ID=FINDS_ALL_READS,Number=1,Type=Integer,Description=\"Detection of SV reads at breakpoint sites, 0 if all reads are found, 1 if not all reads are found\">\n")
                reformat_fo.write("##INFO=<ID=L_INS_LEFT_START,Number=1,Type=Integer,Description=\"Number of reads indicating a long insertion left of the start breakpoint\">\n")
                reformat_fo.write("##INFO=<ID=L_INS_RIGHT_START,Number=1,Type=Integer,Description=\"Number of reads indicating a long insertion right of the start breakpoint\">\n")
                reformat_fo.write("##INFO=<ID=L_INS_LEFT_END,Number=1,Type=Integer,Description=\"Number of reads indicating a long insertion left of the end breakpoint\">\n")
                reformat_fo.write("##INFO=<ID=L_INS_RIGHT_END,Number=1,Type=Integer,Description=\"Number of reads indicating a long insertion right of the end breakpoint\">\n")
                reformat_fo.write("##INFO=<ID=DP_START,Number=1,Type=Integer,Description=\"Coverage at 10 bp interval around SV start position\">\n")
                reformat_fo.write("##INFO=<ID=DP_END,Number=1,Type=Integer,Description=\"Coverage at 10 bp interval around SV end position\">\n")

                reformat_fo.write(line)
            else:
                reformat_fo.write(line)

    vcf_reader = vcf.Reader(open(vcf_file_path, "r"))

    infos_to_add = {"MQ_START":0, "MQ_END":0, "RMDUP_RNAMES":"", "RMDV":0, "MISS_START":0, "MISS_END":0, "STD_START":0, "STD_END":0, "DIST_START":10000, "DIST_END":10000, "FINDS_ALL_READS":0, "L_INS_LEFT_START":0, "L_INS_RIGHT_START":0, "L_INS_LEFT_END":0, "L_INS_RIGHT_END":0, "DP_START":0, "DP_END":0, "MISS_START_NORM":0, "MISS_END_NORM":0}

    for record in vcf_reader:

        filter = "PASS"
        if len(record.FILTER) > 0:
            filter = ",".join(record.FILTER)

        info = []
        for element in record.INFO:
            if type(record.INFO[element]) == list:
                info.append(element + "=" + ",".join([str(x) for x in record.INFO[element]]))
            else:
                info.append(element + "=" + str(record.INFO[element]))

        for info_to_add in infos_to_add:
            if not info_to_add in record.INFO:
                info.append(info_to_add + "=" + str(infos_to_add[info_to_add]))

        correct = True

        try:
            if not "CHR2" in record.INFO:
                info.append("CHR2=" + record.CHROM)
            if not "END" in record.INFO:
                info.append("END=" + str(record.ALT[0].pos))
            if not "SVLEN" in record.INFO:
                info.append("SVLEN=0")
        except AttributeError:
            correct = False

        sample_format = []
        for format in record.FORMAT.split(":"):
            sample_format.append(record.samples[0][format])

        if correct:

            if record.INFO["SVTYPE"] == "BND":
                reformat_fo.write("\t".join(
                    [
                        record.CHROM,
                        str(record.POS),
                        record.ID,
                        record.REF,
                        str(record.ALT[0]),
                        str(record.QUAL),
                        filter,
                        ";".join(info),
                        record.FORMAT,
                        ":".join([str(x) for x in sample_format])
                    ]) + "\n")
            else:
                reformat_fo.write("\t".join(
                    [
                        record.CHROM,
                        str(record.POS),
                        record.ID,
                        record.REF,
                        str(record.ALT[0]),
                        str(record.QUAL),
                        filter,
                        ";".join(info),
                        record.FORMAT,
                        ":".join([str(x) for x in sample_format])
                    ]) + "\n")

def first_stage_filter (vcf_file_path, filtered_vcf_file_path):

    num_variants = 0
    num_passing_variants = 0

    bed_fo = open(filtered_vcf_file_path.replace(".vcf", ".pass.bed"), "w")

    vcf_reader = vcf.Reader(open(vcf_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_vcf_file_path, "w"), template=vcf_reader)
    vcf_pass_writer = vcf.Writer(open(filtered_vcf_file_path.replace(".vcf", ".pass.vcf"), "w"), template=vcf_reader)

    for record in vcf_reader:

        filter = []

        if np.abs(record.INFO["SVLEN"]) < 5000 and record.INFO["SVTYPE"] != "BND":
            filter.append("SVLEN_FILTER")
        if record.samples[0]["DV"] > 15 or record.samples[0]["DV"] < 2:
            filter.append("DV_FILTER")
        if len(record.INFO["RNAMES"]) > 15 or len(record.INFO["RNAMES"]) < 2:
            filter.append("DV_FILTER")

        if record.samples[0]["DR"] > 50:
            filter.append("DR_FILTER")

        """
        is_chr1 = False
        if record.CHROM in ["chr1", "chr2", "chr9", "chr17"]:
            is_chr1 = True
        if record.INFO["CHR2"]== ["chr1", "chr2", "chr9", "chr17"]:
            is_chr1 = True
        if not is_chr1:
            filter.append("chr_FILTER")
        """


        record.FILTER = ",".join(filter)

        num_variants += 1

        if filter == []:
            num_passing_variants += 1
            record.FILTER = "PASS"
            Misc.write_record_as_bed(record, bed_fo)
            vcf_pass_writer.write_record(record)

        vcf_writer.write_record(record)

    vcf_writer.flush()
    vcf_writer.close()

    vcf_pass_writer.flush()
    vcf_pass_writer.close()

    print("Filtering " + str(num_variants) + " variants, " + str(num_passing_variants) + " variants passing")

def second_stage_filter (vcf_file_path, filtered_vcf_file_path):
    num_variants = 0
    num_passing_variants = 0

    bed_fo = open(filtered_vcf_file_path.replace(".vcf", ".pass.bed"), "w")

    vcf_reader = vcf.Reader(open(vcf_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_vcf_file_path, "w"), template=vcf_reader)
    vcf_pass_writer = vcf.Writer(open(filtered_vcf_file_path.replace(".vcf", ".pass.vcf"), "w"), template=vcf_reader)

    for record in vcf_reader:

        filter = []

        if record.INFO["RMDV"] < 2:
            filter.append("RMDV_FILTER")
        if record.INFO["MISS_START"] > 0.11:
            filter.append("MISS_START_FILTER")
        if record.INFO["MISS_END"] > 0.11:
            filter.append("MISS_END_FILTER")
        if record.INFO["MQ_START"] < 20:
            filter.append("MQ_START_FILTER")
        if record.INFO["MQ_END"] < 20:
            filter.append("MQ_END_FILTER")
        if record.INFO["STD_START"] > 100:
            filter.append("STD_START_FILTER")
        if record.INFO["STD_END"] > 100:
            filter.append("STD_END_FILTER")
        if record.INFO["DIST_START"] < 35:
            filter.append("DIST_START_FILTER")
        if record.INFO["DIST_END"] < 35:
            filter.append("DIST_END_FILTER")
        if record.INFO["FINDS_ALL_READS"] == 1:
            filter.append("NOT_ALL_READS_FOUND_FILTER")
        if record.INFO["DP_START"] > 100:
            filter.append("DP_START_FILTER")
        if record.INFO["DP_END"] > 100:
            filter.append("DP_END_FILTER")

        record.FILTER = ",".join(filter)

        num_variants += 1

        if filter == []:
            num_passing_variants += 1
            record.FILTER = "PASS"
            Misc.write_record_as_bed(record, bed_fo)
            vcf_pass_writer.write_record(record)

        vcf_writer.write_record(record)

    vcf_writer.flush()
    vcf_writer.close()

    vcf_pass_writer.flush()
    vcf_pass_writer.close()

    print("Filtering " + str(num_variants) + " variants, " + str(num_passing_variants) + " variants passing")

def third_stage_filter (vcf_file_path, filtered_vcf_file_path):
    num_variants = 0
    num_passing_variants = 0

    bed_fo = open(filtered_vcf_file_path.replace(".vcf", ".pass.bed"), "w")

    vcf_reader = vcf.Reader(open(vcf_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_vcf_file_path, "w"), template=vcf_reader)
    vcf_pass_writer = vcf.Writer(open(filtered_vcf_file_path.replace(".vcf", ".pass.vcf"), "w"), template=vcf_reader)

    for record in vcf_reader:

        filter = []

        if record.INFO["RMDV"] < 2:
            filter.append("GERMLINE_RMDV_FILTER")

        known_germline_variant = filter_known_germline_variants(record, "/home/ubuntu/seq/Miltenyi_ONT/GRCh38_dbVar_common_1000g", 250)
        if known_germline_variant != "":
            filter.append("GERMLINE_DATABASE_FILTER")

        simple_repeat_return_code = filter_for_simple_repeats(record, "/home/ubuntu/data/reference/GRCh38_repeats.sort.simple_repeats.bed")
        if simple_repeat_return_code == 1:
            filter.append("REPEAT_FILTER")

        sniffles_repeat_return_code = filter_for_sniffles_repeats(record, "/home/ubuntu/seq/Miltenyi_ONT/Sniffles/annotations/human_GRCh38_no_alt_analysis_set.trf.pad.cov.frac.bed")
        if sniffles_repeat_return_code == 1:
            filter.append("REPEAT_FILTER")

        if record.INFO["L_INS_LEFT_START"] >= 5 and record.INFO["L_INS_RIGHT_START"] >= 5:
            filter.append("LONG_INSERTION_START_FILTER")
        if record.INFO["L_INS_LEFT_END"] >= 5 and record.INFO["L_INS_RIGHT_END"] >= 5:
            filter.append("LONG_INSERTION_END_FILTER")

        record.FILTER = ",".join(filter)

        num_variants += 1

        if filter == []:
            num_passing_variants += 1
            record.FILTER = "PASS"
            Misc.write_record_as_bed(record, bed_fo)
            vcf_pass_writer.write_record(record)

        vcf_writer.write_record(record)

    vcf_writer.flush()
    vcf_writer.close()

    vcf_pass_writer.flush()
    vcf_pass_writer.close()

    print("Filtering " + str(num_variants) + " variants, " + str(num_passing_variants) + " variants passing")

def exp_filter(vcf_file_path, filtered_vcf_file_path):
    num_variants = 0
    num_passing_variants = 0

    bed_fo = open(filtered_vcf_file_path.replace(".vcf", ".pass.bed"), "w")

    vcf_reader = vcf.Reader(open(vcf_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_vcf_file_path, "w"), template=vcf_reader)
    vcf_pass_writer = vcf.Writer(open(filtered_vcf_file_path.replace(".vcf", ".pass.vcf"), "w"), template=vcf_reader)

    passing_ids = []

    for record in vcf_reader:

        filter = []

        if record.INFO["MISS_START"] / record.INFO["MISS_START_NORM"] > 2.5 and record.INFO["MISS_START"] > 0.05:
            filter.append("MISS_START_FILTER")
        if record.INFO["MISS_END"] / record.INFO["MISS_END_NORM"] > 2.5 and record.INFO["MISS_END"] > 0.05:
            filter.append("MISS_END_FILTER")

        record.FILTER = ",".join(filter)

        num_variants += 1

        if filter == []:
            num_passing_variants += 1
            record.FILTER = "PASS"
            Misc.write_record_as_bed(record, bed_fo)
            vcf_pass_writer.write_record(record)
            passing_ids.append(record.ID)

        vcf_writer.write_record(record)

    vcf_writer.flush()
    vcf_writer.close()

    vcf_pass_writer.flush()
    vcf_pass_writer.close()

    print("Filtering " + str(num_variants) + " variants, " + str(num_passing_variants) + " variants passing")

    cmd_string = "cat "
    for id in passing_ids:
        cmd_string += "assemblies/" + id + ".fasta "
    cmd_string += "> final_passing_assemblies.fasta"
    os.system(cmd_string)

def combine_duplicates (vcf_file_path, filtered_vcf_file_path):
    num_variants = 0
    num_passing_variants = 0

    cutoff = 20

    bed_fo = open(filtered_vcf_file_path.replace(".vcf", ".pass.bed"), "w")

    vcf_reader = vcf.Reader(open(vcf_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_vcf_file_path, "w"), template=vcf_reader)

    vcf_records = []
    for record in vcf_reader:
        vcf_records.append(record)
    copies = np.zeros(len(vcf_records))

    for i in range(0, len(vcf_records)):
        for j in range(0, len(vcf_records)):
            if i != j and copies[i] == 0:
                #Same orientation
                dist_start = 10000
                dist_end = 10000

                if vcf_records[i].CHROM == vcf_records[j].CHROM:
                    dist_start = np.abs(vcf_records[i].POS - vcf_records[j].POS)
                if vcf_records[i].INFO["CHR2"] == vcf_records[j].INFO["CHR2"]:
                    dist_end = np.abs(vcf_records[i].INFO["END"] - vcf_records[j].INFO["END"])

                if dist_start < cutoff and dist_end < cutoff:
                    print(vcf_records[i].ID + " and " + vcf_records[j].ID + " are copies with dist_start " + str(dist_start) + " and dist_end " + str(dist_end))
                    if vcf_records[i].INFO["RMDV"] >= vcf_records[j].INFO["RMDV"]:
                        copies[j] = 1
                    else:
                        copies[i] = 1

                # Different orientation
                dist_start = 10000
                dist_end = 10000

                if vcf_records[i].CHROM == vcf_records[j].INFO["CHR2"]:
                    dist_start = np.abs(vcf_records[i].POS - vcf_records[j].INFO["END"])
                if vcf_records[i].INFO["CHR2"] == vcf_records[j].CHROM:
                    dist_end = np.abs(vcf_records[i].INFO["END"] - vcf_records[j].POS)

                if dist_start < cutoff and dist_end < cutoff:
                    print(vcf_records[i].ID +  " and " + vcf_records[j].ID + " are copies with dist_start " + str(dist_start) + " and dist_end " + str(dist_end))
                    if vcf_records[i].INFO["RMDV"] >= vcf_records[j].INFO["RMDV"]:
                        copies[j] = 1
                    else:
                        copies[i] = 1

    for i in range(0, len(vcf_records)):
        num_variants += 1
        if copies[i] == 0:
            num_passing_variants += 1
            vcf_writer.write_record(vcf_records[i])
            Misc.write_record_as_bed(vcf_records[i], bed_fo)
    vcf_writer.flush()
    vcf_writer.close()

    print("Filtering " + str(num_variants) + " variants, " + str(num_passing_variants) + " variants passing")

def annotate_vcf (vcf_file_path, filtered_vcf_file_path, bam_file_path, logfile_fo):
    logfile_fo.write("Annotating VCF file " + vcf_file_path + "\n")

    vcf_reader = vcf.Reader(open(vcf_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_vcf_file_path, "w"), template=vcf_reader)

    samfile = pysam.AlignmentFile(bam_file_path, "rb")

    for record in vcf_reader:

        print(record.ID)
        logfile_fo.write("\tWorking on variant record " + record.ID + "\n")

        annot_record = analyze_SV_reads(record, samfile, bam_file_path, logfile_fo)
        vcf_writer.write_record(annot_record)

    logfile_fo.flush()
    vcf_writer.flush()
    vcf_writer.close()

def analyze_SV_reads (vcf_record, alignment_file, alignment_file_file_path, logfile_fo):

    # gathering supporting alignments

    finds_all_read_both_breakpoints = False

    logfile_fo.write("\t\tLooking for supporting reads at SV start, expecting " + str(len(vcf_record.INFO["RNAMES"])) + " reads \n")
    supporting_alignments = Read_methods.collect_reads(vcf_record.INFO["RNAMES"], vcf_record.CHROM, vcf_record.POS, 2000, alignment_file, logfile_fo)
    logfile_fo.write("\t\tLooking for supporting reads at SV end, expecting " + str(len(vcf_record.INFO["RNAMES"])) + " reads \n")
    supporting_alignments_end = Read_methods.collect_reads(vcf_record.INFO["RNAMES"], vcf_record.INFO["CHR2"], vcf_record.INFO["END"], 2000, alignment_file, logfile_fo)

    if len(supporting_alignments) == len(vcf_record.INFO["RNAMES"]) and len(supporting_alignments_end) == len(vcf_record.INFO["RNAMES"]):
        finds_all_read_both_breakpoints = True

    if not finds_all_read_both_breakpoints:
        vcf_record.INFO["FINDS_ALL_READS"] = 1

    # Coverage depth
    coverage_bed_file = open("coverage.bed", "w")
    coverage_bed_file.write(vcf_record.CHROM + "\t" + str(vcf_record.POS - 10) + "\t" + str(vcf_record.POS + 10) + "\n")
    coverage_bed_file.write(vcf_record.INFO["CHR2"] + "\t" + str(vcf_record.INFO["END"] - 10) + "\t" + str(vcf_record.INFO["END"] + 10) + "\n")
    coverage_bed_file.close()
    os.system("samtools bedcov coverage.bed " + alignment_file_file_path + " > coverage.st.bed")
    bed_lines = open("coverage.st.bed").readlines()
    vcf_record.INFO["DP_START"] = float(bed_lines[0].split("\t")[3]) / (int(bed_lines[0].split("\t")[2]) - int(bed_lines[0].split("\t")[1]))
    vcf_record.INFO["DP_END"] = float(bed_lines[1].split("\t")[3]) / (int(bed_lines[1].split("\t")[2]) - int(bed_lines[1].split("\t")[1]))

    #Precision

    logfile_fo.write("\t\tLooking for SV precision\n")

    pos_start, std_start = Read_methods.get_end_std_dev(supporting_alignments, logfile_fo)
    pos_end, std_end = Read_methods.get_end_std_dev(supporting_alignments_end, logfile_fo)

    if "PRECISE" in vcf_record.INFO:
        logfile_fo.write("\t\t\tSV is precise\n")
    else:
        logfile_fo.write("\t\t\tSV is imprecise\n")

    logfile_fo.write("\t\t\tStart " + str(pos_start) + "+-" + str(std_start) + "\n")
    logfile_fo.write("\t\t\tEnd " + str(pos_end) + "+-" + str(std_end) + "\n")

    vcf_record.POS = pos_start
    vcf_record.INFO["END"] = pos_end
    if vcf_record.INFO["SVTYPE"] == "BND":
        vcf_record.ALT[0].pos = pos_end
    vcf_record.INFO["STD_START"] = std_start
    vcf_record.INFO["STD_END"] = std_end

    # gathering mean MQ for start and end

    vcf_record.INFO["MQ_START"] = Read_methods.get_reads_mean_MQ(supporting_alignments)
    vcf_record.INFO["MQ_END"] = Read_methods.get_reads_mean_MQ(supporting_alignments_end)

    # gathering mean mismatch fraction for start and end

    vcf_record.INFO["MISS_START"] = Read_methods.get_read_mean_mismatch_fraction(supporting_alignments)
    vcf_record.INFO["MISS_END"] = Read_methods.get_read_mean_mismatch_fraction(supporting_alignments_end)

    norm_miss_start_reads = []
    for norm_read in alignment_file.fetch(vcf_record.CHROM, vcf_record.POS - 100, vcf_record.POS + 100):
        if norm_read.mapping_quality > 0:
            if (norm_read.query_name in vcf_record.INFO["RNAMES"]) == False:
                norm_miss_start_reads.append(norm_read)
    norm_miss_end_reads = []
    for norm_read in alignment_file.fetch(vcf_record.INFO["CHR2"], vcf_record.INFO["END"] - 100, vcf_record.INFO["END"] + 100):
        if norm_read.mapping_quality > 0:
            if (norm_read.query_name in vcf_record.INFO["RNAMES"]) == False:
                norm_miss_end_reads.append(norm_read)

    vcf_record.INFO["MISS_START_NORM"] = Read_methods.get_read_mean_mismatch_fraction(norm_miss_start_reads)
    vcf_record.INFO["MISS_END_NORM"] = Read_methods.get_read_mean_mismatch_fraction(norm_miss_end_reads)


    # writing alignments to file

    fastq_fo = open("fastq/" + vcf_record.ID + ".fastq", "w")
    for r in supporting_alignments:
        Read_methods.write_fastq_original(r, fastq_fo)
    fastq_fo.close()

    logfile_fo.flush()

    # building assembly

    Read_methods.assemble_reads("fastq/" + vcf_record.ID + ".fastq", "/home/ubuntu/seq/Miltenyi_ONT/lamassemble/lamassemble/./lamassemble", "/home/ubuntu/seq/Miltenyi_ONT/lamassemble/lamassemble/train/promethion.mat", vcf_record.ID, "assemblies/" + vcf_record.ID + ".fasta", logfile_fo)

    # mapping reads to assembly and detecting duplicates

    rmdup_rnames, dist_start, dist_end = Read_methods.detect_duplicates(vcf_record, supporting_alignments, "fastq/" + vcf_record.ID + ".fastq", "assemblies/" + vcf_record.ID + ".fasta", "/home/ubuntu/data/tools/minimap2/minimap2", logfile_fo)

    vcf_record.INFO["RMDUP_RNAMES"] = ",".join(rmdup_rnames)
    vcf_record.INFO["RMDV"] = len(rmdup_rnames)

    if dist_start != 10000 and dist_end != 10000:
        vcf_record.INFO["DIST_START"] = dist_start
        vcf_record.INFO["DIST_END"] = dist_end

    logfile_fo.flush()

    return vcf_record

def germline_filter (vcf_file_path, filtered_vcf_file_path, bam_file_path, logfile_fo):

    vcf_reader = vcf.Reader(open(vcf_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_vcf_file_path, "w"), template=vcf_reader)

    samfile = pysam.AlignmentFile(bam_file_path, "rb")

    for record in vcf_reader:
        vcf_writer.write_record(filter_germline_insertions(record, samfile, logfile_fo))

    vcf_writer.flush()
    vcf_writer.close()

def filter_germline_insertions (vcf_record, samfile, logfile_fo):

    logfile_fo.write("Looking at SV " + vcf_record.ID + "\n")

    # looking for read spanning the breakpoint with long (min 200 bp) insertion or by SA

    breakpoints = [[vcf_record.CHROM, vcf_record.POS], [vcf_record.INFO["CHR2"], vcf_record.INFO["END"]]]
    idx = 0

    for breakpoint in breakpoints:
        idx += 1

        logfile_fo.write("\tLooking at breakpoint " + breakpoint[0] + ":" + str(breakpoint[1]) + "\n")

        reads_for_very_long_insertion_left = []
        reads_for_very_long_insertion_right = []

        reads_of_interest = {}
        #for read in samfile.fetch(vcf_record.CHROM, vcf_record.POS - 3000, vcf_record.POS + 3000):
        for read in samfile.fetch(breakpoint[0], breakpoint[1] - 3000, breakpoint[1] + 3000):
            if read.mapping_quality > 0:
                logfile_fo.write("\t\tLooking at read " + read.query_name + "\n")

                # looking for reads spanning the insertion without SA tag
                if vcf_record.POS - read.reference_start > 200 and read.reference_end - vcf_record.POS > 200:

                    has_long_insertion = False

                    cigar_pos = read.reference_start

                    for cigar in read.cigartuples:
                        if cigar[0] in [0, 2]:
                            cigar_pos += cigar[1]
                        if cigar[0] == 1 and cigar[1] > 200:
                            if np.abs(cigar_pos - vcf_record.POS) < 200:
                                has_long_insertion = True
                                logfile_fo.write("\t\t\tRead " + read.query_name + " has an " + str(cigar[1]) + " bp insertion near the SV\n")

                    if has_long_insertion:
                        if read.query_name in reads_of_interest:
                            reads_of_interest[read.query_name].append(read)
                        else:
                            reads_of_interest[read.query_name] = [read]

                # looking for reads with SA tag spanning the deletion
                if read.has_tag("SA"):
                    #logfile_fo.write("\t\tRead " + read.query_name + " has SA tag\n")
                    #logfile_fo.write("\t\tRead " + read.query_name + " has SA tag: " + str(read.get_tag("SA")) + "\n")
                    # at least one SA on the same chromosome
                    SA_on_same_chromosome = False
                    SA_near_breakpoint = False
                    for SA in read.get_tag("SA").split(";"):
                        if SA != "":
                            SA_chrom = SA.split(",")[0]
                            SA_pos = int(SA.split(",")[1])
                            #logfile_fo.write("\t\t\tRead " + read.query_name + " has SA at " + SA_chrom + ":" + str(SA_pos) + "\n")
                            if SA_chrom == read.reference_name:
                                #logfile_fo.write("\t\t\t\tRead " + read.query_name + " has SA on same chromosome with distance " + str(np.abs(SA_pos - breakpoint[1])) + " |" + str(SA_pos) + "-" + str(breakpoint[1]) + "|\n")
                                SA_on_same_chromosome = True
                                #logfile_fo.write("\t\t\t\tDistance of " + read.query_name + " from breakpoint is " + str(np.abs(SA_pos - breakpoint[1])) + " |" + str(SA_pos) + "-" + str(breakpoint[1]) + "|\n")
                                if np.abs(SA_pos - breakpoint[1]) < 20000:
                                    SA_near_breakpoint = True

                            if SA_on_same_chromosome and SA_near_breakpoint:
                                logfile_fo.write("\t\t\tRead " + read.query_name + " has SA near breakpoint\n")
                                if read.query_name in reads_of_interest:
                                    reads_of_interest[read.query_name].append(read)
                                else:
                                    reads_of_interest[read.query_name] = [read]
                            else:
                                p = 0
                                #logfile_fo.write("\t\t\t\tSA is far away from the breakpoints\n")

                # looking for reads with softclipped bases mapping next to the SV
                if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 100:
                    if np.abs(read.reference_start - breakpoint[1]) < 100:
                        logfile_fo.write("\t\t\tRead " + read.query_name + " is right of breakpoint with SA (" + read.reference_name + ":" + str(read.reference_start) + "-" + str(read.reference_end) + ")\n")
                        reads_for_very_long_insertion_right.append(read)
                if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 100:
                    if np.abs(read.reference_end - breakpoint[1]) < 100:
                        logfile_fo.write("\t\t\tRead " + read.query_name + " is left of breakpoint with SA (" + read.reference_name + ":" + str(read.reference_start) + "-" + str(read.reference_end) + ")\n")
                        reads_for_very_long_insertion_left.append(read)


        logfile_fo.write("\t\t*********************\n")
        logfile_fo.write("\t\tFound " + str(len(reads_of_interest)) + " are gathered for further analysis of small or medium insertion\n")
        logfile_fo.write("\t\tFound (" + str(len(reads_for_very_long_insertion_left)) + "/" + str(len(reads_for_very_long_insertion_right)) + ") reads potentially indicating a long insertion at this breakpoint\n")
        logfile_fo.write("\t\t*********************\n")

        if idx == 1:
            vcf_record.INFO["L_INS_LEFT_START"] = len(reads_for_very_long_insertion_left)
            vcf_record.INFO["L_INS_RIGHT_START"] = len(reads_for_very_long_insertion_right)
        if idx == 2:
            vcf_record.INFO["L_INS_LEFT_END"] = len(reads_for_very_long_insertion_left)
            vcf_record.INFO["L_INS_RIGHT_END"] = len(reads_for_very_long_insertion_right)

        # gathering fastq files for insertion assembly and position information for hg38 reference extraction
        germline_insertion_reference_reads_fo = open("germline_fastq/" + vcf_record.ID + ".germline.fastq", "w")
        total_min = vcf_record.POS
        total_max = vcf_record.POS

        num_spanning_reads = 0

        for roi in reads_of_interest:
            start_positions = []
            end_positions = []
            for SA in reads_of_interest[roi]:
                start_positions.append(SA.reference_start)
                end_positions.append(SA.reference_end)
                #bed_out.write(vcf_record.CHROM + "\t" + str(SA.reference_start) + "\t" + str(SA.reference_end) + "\t" + roi + "\n")
            #bed_out.write(vcf_record.CHROM + "\t" + str(np.amin(start_positions)) + "\t" + str(np.amax(end_positions)) + "\t" + roi + "\n")
            if np.amin(start_positions) < total_min:
                total_min = np.amin(start_positions)
            if np.amax(end_positions) > total_max:
                total_max = np.amax(end_positions)

            if breakpoint[1] - np.amin(start_positions) > 100 and np.amax(end_positions) - breakpoint[1] > 100:
                Read_methods.write_fastq_original(reads_of_interest[roi][0], germline_insertion_reference_reads_fo)
                logfile_fo.write("\t\t\tRead " + roi + " spans the SV\n")
                num_spanning_reads += 1

        if num_spanning_reads > 0:

            # Assembly + hg38 extraction + mapping der SV reads auf beide Referenzen gleichzeitig
            germline_insertion_reference_reads_fo.close()
            bed_out = open("germline_intervals/" + vcf_record.ID + ".bed", "w")
            bed_out.write(breakpoint[0] + "\t" + str(total_min) + "\t" + str(total_max) + "\tTotal\n")
            bed_out.flush()

            Read_methods.assemble_reads("germline_fastq/" + vcf_record.ID + ".germline.fastq",
                                        "/home/ubuntu/seq/Miltenyi_ONT/lamassemble/lamassemble/./lamassemble",
                                        "/home/ubuntu/seq/Miltenyi_ONT/lamassemble/lamassemble/train/promethion.mat",
                                        vcf_record.ID + ".germline", "germline_assembly/" + vcf_record.ID + ".germline.fasta", logfile_fo)

            os.system("bedtools getfasta -fi /home/ubuntu/data/reference/Homo_sapiens_assembly38.fasta -bed germline_intervals/" + vcf_record.ID + ".bed > germline_reference/" + vcf_record.ID + ".fasta")
            os.system("cat germline_reference/" + vcf_record.ID + ".fasta germline_assembly/" + vcf_record.ID + ".germline.fasta > combined_reference/" + vcf_record.ID + ".fasta")
            os.system("/home/ubuntu/data/tools/minimap2/minimap2 -ax map-ont -t 5 combined_reference/" + vcf_record.ID + ".fasta fastq/" + vcf_record.ID + ".fastq -Y | samtools view -b - | samtools sort - > germline_mappings/" + vcf_record.ID + ".bam")
            os.system("samtools index germline_mappings/" + vcf_record.ID + ".bam")

            #Schauen und zaehlen ob die read auf hg38 oder auf die assembly mappen

            comparative_samfile = pysam.AlignmentFile("germline_mappings/" + vcf_record.ID + ".bam", "rb")

            rmdv = vcf_record.INFO["RMDV"]
            rmdv_rnames = vcf_record.INFO["RMDUP_RNAMES"]

            for read in comparative_samfile:
                if not read.is_secondary and not read.is_supplementary and read.is_mapped:
                    logfile_fo.write("\t\tRead " + read.query_name + " maps on " + read.reference_name + "\n")
                    if read.reference_name == vcf_record.ID + ".germline":
                        try:
                            rmdv_rnames.remove(read.query_name)
                            logfile_fo.write("\t\t\tDeleting read " + read.query_name + " from RMDV, " + str(len(rmdv_rnames)) + " reads supporting the SV remaining\n")
                        except ValueError:
                            logfile_fo.write("\t\t\t" + read.query_name + " not found in RMDV list\n")

            if len(rmdv_rnames) < vcf_record.INFO["RMDV"]:
                vcf_record.INFO["RMDV"] = len(rmdv_rnames)
                vcf_record.INFO["RMDUP_RNAMES"] = rmdv_rnames

        else:
            logfile_fo.write("\t\t*********************\n")
            logfile_fo.write("\t\tNo germline insertion at " + breakpoint[0] + ":" + str(breakpoint[1]) + "\n")
            logfile_fo.write("\t\t*********************\n")

    return vcf_record

def filter_known_germline_variants (vcf_record, dbVar_common_file_file_path, tolerance):
    #print("\tLooking for known germline SVs explaining " + vcf_record.ID + "\n")
    if vcf_record.INFO["SVTYPE"] != "BND":
        for line in open(dbVar_common_file_file_path):
            if not line.startswith("#chrom"):
                chrom = line.split("\t")[0]
                start = int(line.split("\t")[1])
                end = int(line.split("\t")[2])

                if vcf_record.CHROM == chrom:
                    start_diff = np.abs(vcf_record.POS - start)
                    end_diff = np.abs(vcf_record.INFO["END"] - end)
                    if (start_diff + end_diff) < tolerance:
                        print("\tFound germline SV " + line.split("\t")[3] + " explaining this SV " + vcf_record.ID)
                        return line.split("\t")[3]
        return ""
    return ""

def filter_for_simple_repeats (vcf_record, repeats_file_path):
    i = 0
    start_distance_to_next_repeat = 10000
    end_distance_to_next_repeat = 10000
    for breakpoint in [[vcf_record.CHROM, vcf_record.POS], [vcf_record.INFO["CHR2"], vcf_record.INFO["END"]]]:
        i += 1
        min_dist = 100000
        for line in open(repeats_file_path):
            if line.split("\t")[3].startswith("("):
                if line.split("\t")[0] == breakpoint[0]:
                    # print(line.rstrip())
                    dist = np.amin([np.abs(breakpoint[1] - int(line.split("\t")[1])),
                                    np.abs(breakpoint[1] - int(line.split("\t")[2]))])
                    # print(dist)
                    if breakpoint[1] >= int(line.split("\t")[1]) and breakpoint[1] <= int(line.split("\t")[2]):
                        dist = 0
                    if dist < min_dist:
                        min_dist = dist

        if i == 1:
            start_distance_to_next_repeat = min_dist
        if i == 2:
            end_distance_to_next_repeat = min_dist

    if end_distance_to_next_repeat < 20 and start_distance_to_next_repeat < 20:
        print("\tBoth breakpoints of SV " + vcf_record.ID + " are in simple repeats")
        return 1
    else:
        return 0

def filter_for_sniffles_repeats (vcf_record, repeats_file_path):
    i = 0
    start_distance_to_next_repeat = 10000
    end_distance_to_next_repeat = 10000
    for breakpoint in [[vcf_record.CHROM, vcf_record.POS], [vcf_record.INFO["CHR2"], vcf_record.INFO["END"]]]:
        i += 1
        min_dist = 100000
        for line in open(repeats_file_path):
            if float(line.split("\t")[3]) > 1.8:
                if line.split("\t")[0] == breakpoint[0]:
                    # print(line.rstrip())
                    dist = np.amin([np.abs(breakpoint[1] - int(line.split("\t")[1])),
                                    np.abs(breakpoint[1] - int(line.split("\t")[2]))])
                    # print(dist)
                    if breakpoint[1] >= int(line.split("\t")[1]) and breakpoint[1] <= int(line.split("\t")[2]):
                        dist = 0
                    if dist < min_dist:
                        min_dist = dist

        if i == 1:
            start_distance_to_next_repeat = min_dist
        if i == 2:
            end_distance_to_next_repeat = min_dist

    if end_distance_to_next_repeat < 20 and start_distance_to_next_repeat < 20:
        print("\tBoth breakpoints of SV " + vcf_record.ID + " are in sniffles repeats")
        return 1
    else:
        return 0
