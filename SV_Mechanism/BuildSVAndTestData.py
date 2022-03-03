import General.LoadData
import General.ChromosomeMetrics

def build_BED_around_SV_breakpoints (sv_bed_file_path):
    sv_bed_data = General.LoadData.load_sv_bed(sv_bed_file_path)

    print(sv_bed_data)