import pandas as pd

# Loads data from bed file with all SV informations and output breakpoints per SV in a dictionary with the SV id as key

def load_sv_bed (sv_bed_file_path):

    sv_dict = {}

    for line in open(sv_bed_file_path):
        chrom = line.rstrip().split("\t")[0]
        start = int(line.rstrip().split("\t")[1])
        end = int(line.rstrip().split("\t")[2])
        descr = line.rstrip().split("\t")[3]
        id = descr.split("_")[0]

        svtype = "NULL"
        if ":" in descr:
            svtype = "BND"

            if id in sv_dict:
                sv_dict[id].append([chrom, start, svtype])
            else:
                sv_dict[id] = [[chrom, start, svtype]]

        else:
            svtype = descr.split("_")[1]

            if id in sv_dict:
                sv_dict[id].append([chrom, start, svtype])
                sv_dict[id].append([chrom, end, svtype])
            else:
                sv_dict[id] = [[chrom, start, svtype]]
                sv_dict[id].append([chrom, end, svtype])
    """
    # Potential output
    for sv in sv_dict:
        print(sv)
        for item in sv_dict[sv]:
            print(item)
        print("-----------------")
    """
    
    return sv_dict

# Loads data from bed file with all SV informations and output SV id and SV type

def load_sv_bed_agg_by_id (sv_bed_file_path):
    sv_dict = load_sv_bed(sv_bed_file_path)
    agg_sv_dict = {}
    for sv in sv_dict:
        agg_sv_dict[sv] = sv_dict[sv][0][-1]
    return agg_sv_dict

# reads fasta file with one or multiple sequences, returns array mit [[seq, name], ...]

def readFASTA (fasta_file_path):
    seqs = []

    fasta_name = ""
    fasta_seq = ""
    for line in open(fasta_file_path):
        if line.startswith(">"):
            if fasta_name == "":
                fasta_name = line.rstrip().replace(">", "")
            else:
                seqs.append([fasta_name, fasta_seq])
                fasta_name = line.rstrip().replace(">", "")
                fasta_seq = ""
        else:
            fasta_seq += line.rstrip()
    if fasta_name != "":
        seqs.append([fasta_name, fasta_seq])

    return seqs

# load data from GATK4 CNV seg-File as [CHROM, START, END, N_COV_DATAPOINTS, N_HET_SNPS, CNR_10, CNR, CNR_90, AF_10, AF, AF_90]

def load_seg_data_GATK_CNV_seg_file(seg_file):
    seg_data = pd.DataFrame()

    for line in open(seg_file):

        if line.startswith("@") == False and line.startswith("CONTIG") == False:
            chrom = line.split("\t")[0]
            start = int(line.split("\t")[1])
            end = int(line.split("\t")[2])
            n_cnr_datapoints = int(line.split("\t")[3])
            n_hets_snps = int(line.split("\t")[4])

            cnr_10 = float(line.split("\t")[5])
            cnr_50 = float(line.split("\t")[6])
            cnr_90 = float(line.split("\t")[7])

            af_10 = float(line.split("\t")[8])
            af_50 = float(line.split("\t")[9])
            af_90 = float(line.split("\t")[10])

            seg_data = seg_data.append(
                pd.DataFrame({"CHROM": [chrom],
                              "START": [start],
                              "END": [end],
                              "N_COV_DATAPOINTS": [n_cnr_datapoints],
                              "N_HET_SNPS": [n_hets_snps],
                              "CNR_10": [cnr_10],
                              "CNR": [cnr_50],
                              "CNR_90": [cnr_90],
                              "AF_10": [af_10],
                              "AF": [af_50],
                              "AF_90": [af_90]
                              }), ignore_index=True
            )
    return seg_data

# load data from GATK4 SNV TSV file as [CHROM, POS, DP, AF], build GATK4 SNV TSV with gatk VariantsToTable -F CHROM -F POS -F FILTER -GF AF -GF AD

def load_SNV_data_GATK_SNV_table_file(table_file, tumor_id):
    snv_data = pd.read_csv(table_file, sep="\t")
    return_snv_data = pd.DataFrame()

    chrom_array = []
    pos_array = []
    dp_array = []
    af_array = []

    tumor_af_index = -1
    tumor_ad_index = -1

    for i in range(0, len(snv_data.columns)):
        if snv_data.columns[i] == tumor_id + ".AF":
            tumor_af_index = i
        if snv_data.columns[i] == tumor_id + ".AD":
            tumor_ad_index = i

    for line in open(table_file):

        if line.startswith("CHROM") == False:

            try:

                chrom = line.split("\t")[0]
                pos = int(line.split("\t")[1])

                ad = line.split("\t")[tumor_ad_index].split(",")
                dp = int(ad[0]) + int(ad[1])
                af = float(line.split("\t")[tumor_af_index])

                chrom_array.append(chrom)
                pos_array.append(pos)
                dp_array.append(dp)
                af_array.append(af)

            except ValueError:
                ve = 0

    return_snv_data = return_snv_data.append(
        pd.DataFrame({"CHROM": chrom_array,
                      "POS": pos_array,
                      "DP": dp_array,
                      "AF": af_array
                      }), ignore_index=True
    )
    return return_snv_data
