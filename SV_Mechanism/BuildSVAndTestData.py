import General.LoadData
import General.ChromosomeMetrics

# Writes BED file with padding around SV breakpoints, inclusion of templates from templated insertion is optional

def build_BED_around_SV_breakpoints (sv_bed_file_path, padded_bed_file_path, include_templates, padding):
    sv_bed_data = General.LoadData.load_sv_bed(sv_bed_file_path)
    fo_padded_bed_file = open(padded_bed_file_path, "w")

    for sv_id in sv_bed_data:
        for breakpoint in sv_bed_data[sv_id]:
            include_sv = True
            if (not include_templates) and breakpoint[2].startswith("template"):
                include_sv = False

            if include_templates:
                fo_padded_bed_file.write("\t".join(
                    [
                        breakpoint[0],
                        str(breakpoint[1] - padding),
                        str(breakpoint[1] + padding),
                        sv_id
                    ]
                ) + "\n")

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