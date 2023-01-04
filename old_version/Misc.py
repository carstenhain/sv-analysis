import vcf
import os
import Read_methods

def write_record_as_bed (record, bed_fo):

    if record.INFO["SVTYPE"] =="BND":
        bed_fo.write("\t".join([record.CHROM, str(record.POS), str(record.POS + 1), str(record.ID) + "_" + record.INFO["CHR2"] + ":" + str(record.ALT[0].pos)]) + "\n")
        bed_fo.write("\t".join([record.INFO["CHR2"], str(record.ALT[0].pos), str(record.ALT[0].pos + 1), str(record.ID) + "_" + record.CHROM + ":" + str(record.POS)]) + "\n")
    else:
        bed_fo.write("\t".join([record.CHROM, str(record.POS), str(record.INFO["END"]), str(record.ID) + "_" + record.INFO["SVTYPE"]]) + "\n")


def get_repeat_fasta (repeat_name, min_length, coloring_cutoff, num, fasta_path):
    num_repeats_saved = 0
    tmp_bed = open("tmp.bed", "w")
    for line in open("/home/ubuntu/data/reference/GRCh38_repeats.sort.bed"):
        length = int(line.split("\t")[2]) - int(line.split("\t")[1])
        name = line.split("\t")[3]
        coloring = int(line.split("\t")[8].split(",")[0])
        if length > min_length and coloring < coloring_cutoff and num_repeats_saved < num:
            tmp_bed.write("\t".join(line.split("\t")[0:3]) + "\n")
            num_repeats_saved += 1
    tmp_bed.close()
    os.system("bedtools getfasta -fi /home/ubuntu/data/reference/Homo_sapiens_assembly38.fasta -bed tmp.bed > " + fasta_path)
    os.system("rm tmp.bed")
    tmp_log = open("tmp.log.txt", "w")
    Read_methods.assemble_reads(fasta_path,
                                "/home/ubuntu/seq/Miltenyi_ONT/lamassemble/lamassemble/./lamassemble",
                                "/home/ubuntu/seq/Miltenyi_ONT/lamassemble/lamassemble/train/promethion.mat",
                                fasta_path.replace(".fasta", ""), fasta_path.replace(".fasta", ".assembly.fasta"), tmp_log)
    tmp_log.close()
    os.system("rm tmp.log.txt")
