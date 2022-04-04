import BuildSVAndTestData
import General.LoadData
import General.SequenceMethods
import numpy as np

# Padding 25 for RAG heptamer analysis
"""
BuildSVAndTestData.build_BED_around_SV_breakpoints("C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WGS\\SV\\CD3Iso.Manuell.Truth_Set.plus.Add_In.templates.manuell_kontrolliert.Final.Sortiert.bed", "..\\Resources\\SV.breakpoints.padding.25.bed", True, 25)
BuildSVAndTestData.build_BED_around_random_breakpoints(200, "..\\Resources\\Random.breakpoints.padding.25.bed", 25)
"""
# Padding 12 for microhomology analysis
"""
BuildSVAndTestData.build_BED_around_SV_breakpoints("C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WGS\\SV\\CD3Iso.Manuell.Truth_Set.plus.Add_In.templates.manuell_kontrolliert.Final.Sortiert.bed", "..\\Resources\\SV.breakpoints.padding.12.bed", True, 12)
BuildSVAndTestData.build_BED_around_random_breakpoints(200, "..\\Resources\\Random.breakpoints.padding.12.bed", 12)
"""
# bedtools getfasta -fi hg38.fasta -bed padded.bed -fo padded.fasta
# FIMO for the padding 25 fasta files, p-value < 1

# Microhomology data for the SVs

data_file = open("Data\\SV_microhomology.txt", "w")
SVs = General.LoadData.load_sv_bed_agg_by_id("C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WGS\\SV\\CD3Iso.Manuell.Truth_Set.plus.Add_In.templates.manuell_kontrolliert.Final.Sortiert.bed")
fasta_data = BuildSVAndTestData.aggregate_FASTA("..\\Resources\\SV.breakpoints.padding.12.fasta")
for sv_id in fasta_data:
    if len(fasta_data[sv_id]) == 2:

        r_normal = [General.SequenceMethods.getMicrohomologyBaseCount(fasta_data[sv_id][0], fasta_data[sv_id][1]),
                    "normal"]
        r_complement = [General.SequenceMethods.getMicrohomologyBaseCount(fasta_data[sv_id][0],
                                                                          General.SequenceMethods.getComplement(
                                                                              fasta_data[sv_id][1])), "complement"]
        r_reverse = [General.SequenceMethods.getMicrohomologyBaseCount(fasta_data[sv_id][0],
                                                                       General.SequenceMethods.getReverse(
                                                                           fasta_data[sv_id][1])), "reverse"]
        r_reversecomplement = [General.SequenceMethods.getMicrohomologyBaseCount(fasta_data[sv_id][0],
                                                                                 General.SequenceMethods.getReverseComplement(
                                                                                     fasta_data[sv_id][1])),
                               "reverse complement"]

        max_microhomology_strech = 0
        if SVs[sv_id] == "BND":
            max_microhomology_strech = np.amax([r_normal[0][4], r_complement[0][4], r_reverse[0][4], r_reversecomplement[0][4]])
        if SVs[sv_id] == "DEL" or SVs[sv_id] == "DUP":
            max_microhomology_strech = r_normal[0][4]
        if SVs[sv_id] == "INV":
            max_microhomology_strech = r_reversecomplement[0][4]

        data_file.write(str(max_microhomology_strech) + "\t" + str(sv_id) + "\n")
data_file.close()

data_file = open("Data\\Random_microhomology.txt", "w")
fasta_data = BuildSVAndTestData.aggregate_FASTA("..\\Resources\\Random.breakpoints.padding.12.fasta")
ids = []
for sv_id in fasta_data:
    ids.append(sv_id)
for i in range(0, int(len(ids) / 2)):
    r_normal = [General.SequenceMethods.getMicrohomologyBaseCount(fasta_data[ids[2 * i]][0], fasta_data[ids[2 * i + 1]][0]), "normal"]
    r_complement = [General.SequenceMethods.getMicrohomologyBaseCount(fasta_data[ids[2 * i]][0],
                                                                      General.SequenceMethods.getComplement(fasta_data[ids[2 * i + 1]][0])), "complement"]
    r_reverse = [General.SequenceMethods.getMicrohomologyBaseCount(fasta_data[ids[2 * i]][0],
                                                                   General.SequenceMethods.getReverse(fasta_data[ids[2 * i + 1]][0])), "reverse"]
    r_reversecomplement = [General.SequenceMethods.getMicrohomologyBaseCount(fasta_data[ids[2 * i]][0],
                                                                             General.SequenceMethods.getReverseComplement(fasta_data[ids[2 * i + 1]][0])), "reverse complement"]

    max_microhomology_strech = np.amax([r_normal[0][4], r_complement[0][4], r_reverse[0][4]])

    data_file.write(str(max_microhomology_strech) + "\n")
data_file.close()


# LCR distance

sv_dict = General.LoadData.load_sv_bed("C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WGS\\SV\\CD3Iso.Manuell.Truth_Set.plus.Add_In.templates.manuell_kontrolliert.Final.Sortiert.bed")
"""
data_file = open("Data\\SV_LCR_distance.txt", "w")
for sv in sv_dict:
    print(sv)
    for item in sv_dict[sv]:
        data_file.write("\t".join([str(x) for x in item]))
        lcr_dist, lcr_name, lcr_dir = General.SequenceMethods.getDistanceToNextLCR(item[0], item[1], "D:\\GRCh38_repeats.sort.bed")
        data_file.write("\t" + sv + "\t" + str(lcr_dist) + "\t" + lcr_name + "\t" + lcr_dir + "\n")
        data_file.flush()

"""
sv_dict = {}
i = 0
for line in open("..\\Resources\\Random.breakpoints.padding.12.bed"):
    id = str(int(i / 2))
    if id in sv_dict:
        sv_dict[id].append([line.split("\t")[0], int(line.split("\t")[1]) + 12, "SV"])
    else:
        sv_dict[id] = [[line.split("\t")[0], int(line.split("\t")[1]) + 12, "SV"]]
    i += 1

data_file = open("Data\\Random_LCR_distance.txt", "w")
print(len(sv_dict))
for sv in sv_dict:
    print(sv)
    for item in sv_dict[sv]:
        data_file.write("\t".join([str(x) for x in item]))
        lcr_dist, lcr_name, lcr_dir = General.SequenceMethods.getDistanceToNextLCR(item[0], item[1], "D:\\GRCh38_repeats.sort.bed")
        data_file.write("\t" + sv + "\t" + str(lcr_dist) + "\t" + lcr_name + "\t" + lcr_dir + "\n")
        data_file.flush()