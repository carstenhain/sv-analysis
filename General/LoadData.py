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

# # Loads data from bed file with all SV informations and output SV id and SV type

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

