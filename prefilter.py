import vcf
import pysam
import numpy as np

def get_mean_MQ (reads):
    data = []
    for read in reads:
        data.append(read.mapping_quality)
    return np.mean(data)

def get_mean_pID (reads):
    data = []
    for read in reads:
        data.append((1 - read.get_tag("NM") / (read.reference_end - read.reference_start)) * 100)
    return np.mean(data)

# annotates SV record with MQ and PID values for both breakpoints
# record : vcf record of this sv
# samfile : pysam samfile of the sample
# padding: search area around the breakpoint for finding SV supporting reads
def annotate_record (record, samfile, padding):


    if record.INFO["SVTYPE"] == "BND":
        breakpoints = \
            [
                [record.CHROM, record.POS, "START"],
                [record.ALT[0].chr, record.ALT[0].pos, "END"]
            ]
    else:
        breakpoints = \
            [
                [record.CHROM, record.POS, "START"],
                [record.INFO["CHR2"], record.INFO["END"], "END"]
            ]

    for b in breakpoints:

        sv_reads = []
        normal_reads = []

        for read in samfile.fetch(b[0], np.amax([0, b[1] - padding]), b[1] + padding):
            if read.query_name in record.INFO["RNAMES"]:
                sv_reads.append(read)
            else:
                normal_reads.append(read)

        record.INFO["MQ_SV_" + b[2]] = get_mean_MQ(sv_reads)
        record.INFO["MQ_NORM_" + b[2]] = get_mean_MQ(normal_reads)
        record.INFO["PID_SV_" + b[2]] = get_mean_pID(sv_reads)
        record.INFO["PID_NORM_" + b[2]] = get_mean_pID(normal_reads)

    return record

# write FILTER fields depending on MQ (<20) and PID (diff > 7) values
# record : vcf record of this sv
def filter_record (record):

    # remove all previous filters
    record.FILTER = []

    # filter for MQ < 20 at start and/or end
    if np.any([record.INFO["MQ_SV_START"] < 20, record.INFO["MQ_SV_END"] < 20]):
        record.FILTER.append("MQ")

    # filter for PID
    if np.any([record.INFO["PID_SV_START"] + 7 < record.INFO["PID_NORM_START"], record.INFO["PID_SV_END"] + 7 < record.INFO["PID_NORM_END"]]):
        record.FILTER.append("PID")

    return record

# adds info and filter lines to the vcf header
def modify_vcf_header (vcf_file_path, info_vcf_file_path):
    fo = open(info_vcf_file_path, "w")
    vcf_file_lines = open(vcf_file_path).readlines()

    for idx, line in enumerate(vcf_file_lines):
        if line.startswith("#"):
            if idx > 1:
                # add info fields
                if vcf_file_lines[idx - 1].startswith("##INFO") and not vcf_file_lines[idx].startswith("##INFO"):
                    # add info fields for mq
                    fo.write("##INFO=<ID=MQ_SV_START,Number=1,Type=Float,Description=\"MQ of SV supporting reads at start breakpoint\">\n")
                    fo.write("##INFO=<ID=MQ_SV_END,Number=1,Type=Float,Description=\"MQ of SV supporting reads at end breakpoint\">\n")
                    fo.write("##INFO=<ID=MQ_NORM_START,Number=1,Type=Float,Description=\"MQ of normal reads at start breakpoint\">\n")
                    fo.write("##INFO=<ID=MQ_NORM_END,Number=1,Type=Float,Description=\"MQ of normal reads at end breakpoint\">\n")

                    # add info fields for pid
                    fo.write("##INFO=<ID=PID_SV_START,Number=1,Type=Float,Description=\"Percent identity of SV supporting reads at start breakpoint\">\n")
                    fo.write("##INFO=<ID=PID_SV_END,Number=1,Type=Float,Description=\"Percent identity of SV supporting reads at end breakpoint\">\n")
                    fo.write("##INFO=<ID=PID_NORM_START,Number=1,Type=Float,Description=\"Percent identity of normal reads at start breakpoint\">\n")
                    fo.write("##INFO=<ID=PID_NORM_END,Number=1,Type=Float,Description=\"Percent identity of normal reads at start breakpoint\">\n")

                # add filters
                if vcf_file_lines[idx - 1].startswith("##FILTER") and not vcf_file_lines[idx].startswith("##FILTER"):
                    fo.write("##FILTER=<ID=MQ,Description=\"SV reads have lower MQ than normal reads\">\n")
                    fo.write("##FILTER=<ID=PID,Description=\"SV reads have lower percent identity than normal reads\">\n")


        fo.write(line)

# applies modify_vcf_header and annotate_record to a complete vcf file
def prefilter (vcf_file_path, sam_file_path, data_folder):

    total_variants = 0
    short_variants = 0
    filtered_variants = 0
    passing_variants = 0

    # add info fields and filter to vcf header
    modify_vcf_header(vcf_file_path, data_folder + vcf_file_path.split(".vcf")[0].split("/")[-1] + ".annot.vcf")

    # load vcf files
    vcf_reader = vcf.Reader(open(data_folder + vcf_file_path.split(".vcf")[0].split("/")[-1] + ".annot.vcf", "r"))
    vcf_writer = vcf.Writer(open(data_folder + vcf_file_path.split(".vcf")[0].split("/")[-1] + ".prefilter.vcf", "w"), template=vcf_reader)

    # load matching samfile
    samfile = pysam.AlignmentFile(sam_file_path, "rb")

    # add prefilter annotations and filter record
    for record in vcf_reader:

        total_variants += 1

        """
        temporary fix of filtering SVs smaller than 5000 bp
        """
        remove_sv = True
        if record.INFO["SVTYPE"] == "BND":
            remove_sv = False
        if np.abs(record.INFO["SVLEN"]) >= 5000:
            remove_sv = False

        if not remove_sv:
            prefiltered_record = filter_record(annotate_record(record, samfile, 2000))
            vcf_writer.write_record(prefiltered_record)
            if prefiltered_record.FILTER == []:
                passing_variants += 1
            else:
                filtered_variants += 1
        else:
            short_variants += 1

    vcf_writer.close()

    return total_variants, short_variants, filtered_variants, passing_variants
