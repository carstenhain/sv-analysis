import General.LoadData

def build_absolute_files (tumor_id, seg_file_path, snv_tsv_file_path, absolute_cnv_file_path, absolute_maf_file_path):

    # CNV data file

    cnv_data_tumor = General.LoadData.load_seg_data_GATK_CNV_seg_file(seg_file_path)
    cnv_data_tumor_file = open(absolute_cnv_file_path, "w")
    cnv_data_tumor_file.write("\t".join(["Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean"]) + "\n")

    for i in range(0, len(cnv_data_tumor["CHROM"])):
        chrom = cnv_data_tumor["CHROM"][i]
        if chrom != "chrX" and chrom != "chrY":
            cnv_data_tumor_file.write("\t".join(
                [
                    tumor_id,
                    cnv_data_tumor["CHROM"][i].replace("chr", ""),
                    str(cnv_data_tumor["START"][i]),
                    str(cnv_data_tumor["END"][i]),
                    str(cnv_data_tumor["N_COV_DATAPOINTS"][i]),
                    str(cnv_data_tumor["CNR"][i])
                ]
            ) + "\n")

    cnv_data_tumor_file.close()

    # SNV data file

    snv_data_tumor_file = open(absolute_maf_file_path, "w")
    snv_data_tumor_file.write("\t".join(
        ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_position", "End_position",
         "Strand", "dbSNP_Val_Status", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "t_alt_count",
         "t_ref_count"]) + "\n")
    snv_data_tumor = General.LoadData.load_SNV_data_GATK_SNV_table_file(snv_tsv_file_path, tumor_id)

    for i in range(0, len(snv_data_tumor["CHROM"])):
        chrom = snv_data_tumor["CHROM"][i]
        if chrom != "chrX" and chrom != "chrY":
            dp = snv_data_tumor["DP"][i]
            af = snv_data_tumor["AF"][i]
            tum = int(dp * af)
            ref = dp - tum

            snv_data_tumor_file.write("\t".join(
                [
                    ".",
                    ".",
                    "project",
                    "38",
                    snv_data_tumor["CHROM"][i].replace("chr", ""),
                    str(snv_data_tumor["POS"][i]),
                    str(snv_data_tumor["POS"][i]),
                    "+",
                    "",
                    tumor_id,
                    ".",
                    str(tum),
                    str(ref)
                ]
            ) + "\n")

    snv_data_tumor_file.close()

build_absolute_files("CD3Iso", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\SEG\\CD3Iso_final.modelFinal.seg", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\VCF\\CD3Iso_final-filtered.PASS.snv.tsv", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\ABSOLUTE\\CD3Iso.ABSOLUTE.CNV.txt", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\ABSOLUTE\\CD3Iso.ABSOLUTE.SNV.txt")
build_absolute_files("CD4Iso", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\SEG\\CD4Iso_final.modelFinal.seg", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\VCF\\CD4Iso_final-filtered.PASS.snv.tsv", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\ABSOLUTE\\CD4Iso.ABSOLUTE.CNV.txt", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\ABSOLUTE\\CD4Iso.ABSOLUTE.SNV.txt")
build_absolute_files("CD3Durch", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\SEG\\CD3Durch_final.modelFinal.seg", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\VCF\\CD3Durch_final-filtered.PASS.snv.tsv", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\ABSOLUTE\\CD3Durch.ABSOLUTE.CNV.txt", "C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\ABSOLUTE\\CD3Durch.ABSOLUTE.SNV.txt")

r_script = open("C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WXS\\WXS.absolute.rscript.txt", "w")
r_script.write("library(DoAbsolute)\nexample_path = system.file(\"extdata\", package = \"DoAbsolute\", mustWork = T)\nlibrary(data.table)\n")

segs = []
mafs = []

for tumor_id in ["CD3Iso", "CD4Iso", "CD3Durch"]:

    r_script.write("\n".join(
        [
            "seg_" + tumor_id + " = file.path(example_path, \"" + tumor_id + ".ABSOLUTE.CNV.txt\")",
            "maf_" + tumor_id + " = file.path(example_path, \"" + tumor_id + ".ABSOLUTE.SNV.txt\")",
            "seg_" + tumor_id + " = fread(seg_" + tumor_id + ")",
            "maf_" + tumor_id + " = fread(maf_" + tumor_id + ")"
        ]
    ) + "\n")

    segs.append("seg_" + tumor_id)
    mafs.append("maf_" + tumor_id)

r_script.write("\nSeg = Reduce(rbind, list(" + ",".join(segs) + "))\n")
r_script.write("Maf = Reduce(rbind, list(" + ",".join(mafs) + "))\n")

r_script.write("Seg$Sample = substr(Seg$Sample, 1, 15)\n")
r_script.write("Maf$Tumor_Sample_Barcode = substr(Maf$Tumor_Sample_Barcode, 1, 15)\n")
r_script.write("DoAbsolute(Seg = Seg, Maf = Maf, platform = \"Illumina_WES\", copy.num.type = \"total\", results.dir = \"phs000725\", nThread = 6, keepAllResult = TRUE, verbose = TRUE)\n")

r_script.close()