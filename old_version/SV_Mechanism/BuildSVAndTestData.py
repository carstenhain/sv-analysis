import General.LoadData
import General.ChromosomeMetrics

# Writes BED file with padding around SV breakpoints, inclusion of templates from templated insertion is optional

def build_BED_around_SV_breakpoints (sv_bed_file_path, padded_bed_file_path, include_templates, padding):
    sv_bed_data = General.LoadData.load_sv_bed(sv_bed_file_path)
    fo_padded_bed_file = open(padded_bed_file_path, "w")

    for sv_id in sv_bed_data:
        i = 0
        for breakpoint in sv_bed_data[sv_id]:
            include_sv = True
            if (not include_templates) and breakpoint[2].startswith("template"):
                include_sv = False

            if include_sv:
                fo_padded_bed_file.write("\t".join(
                    [
                        breakpoint[0],
                        str(breakpoint[1] - padding),
                        str(breakpoint[1] + padding),
                        sv_id + "_" + str(i)
                    ]
                ) + "\n")
                i += 1

    fo_padded_bed_file.flush()

# Writes BED file with padding around random genomic locations

def build_BED_around_random_breakpoints (n_points, padded_bed_file_path, padding):
    fo_padded_bed_file = open(padded_bed_file_path, "w")

    for i in range(0, n_points):
        chrom, pos = General.ChromosomeMetrics.get_random_position("..\\Resources\\karyotype.human.hg38.txt")
        fo_padded_bed_file.write("\t".join(
            [
                chrom,
                str(pos - padding),
                str(pos + padding)
            ]
        ) + "\n")

    fo_padded_bed_file.flush()

# Reads FIMO TSV output saves highest score for each breakpoint

def read_FIMO_output (fimo_output_file_path):
    data = {}

    for line in open(fimo_output_file_path):

        if line.startswith("1"):

            seq_name = line.split("\t")[2]
            score = float(line.split("\t")[6])

            if seq_name in data:
                if score > data[seq_name]:
                    data[seq_name] = score
            else:
                data[seq_name] = score

    data_per_sv = {}

    for breakpoint in data:
        sv_id = breakpoint.split("_")[0]
        if sv_id in data_per_sv:
            data_per_sv[sv_id].append(data[breakpoint])
        else:
            data_per_sv[sv_id] = [data[breakpoint]]

    return data_per_sv

# Aggregates multiple FASTA by SV ID

def aggregate_FASTA (fasta_file_path):
    fasta_data = {}
    raw_fasta_data = General.LoadData.readFASTA(fasta_file_path)

    for seq in raw_fasta_data:
        sv_id = seq[0].split("_")[0]
        breakpoint_seq = seq[1]
        if sv_id in fasta_data:
            fasta_data[sv_id].append(breakpoint_seq)
        else:
            fasta_data[sv_id] = [breakpoint_seq]

    return fasta_data

# Reads SVType from
def getSVIDs (sv_bed_file_path):
    ids = []
    svtypes = []

    for line in open(sv_bed_file_path):
        id = line.split("\t")[3].split("_")[0]
        svtype = "-1"
        if "DEL" in line.split("\t")[3]:
            svtype = "DEL"
        if "INV" in line.split("\t")[3]:
            svtype = "INV"
        if "DUP" in line.split("\t")[3]:
            svtype = "DUP"
        if ":" in line.split("\t")[3]:
            svtype = "BND"
        if not (id in ids):
            ids.append(id)
            svtypes.append(svtype)

    return ids, svtypes
