import numpy as np

def getComplement (seq):
    complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
    new_seq = ""
    for i in range(0, len(seq)):
        new_seq += complement[seq[i]]
    return new_seq

def getReverse (seq):
    reverse = ""
    for i in range(0, len(seq)):
        reverse += seq[len(seq) - 1 - i]
    return reverse

def getReverseComplement(seq):
    return getReverse(getComplement(seq))

# Calculates overlap between sequences a and b (len (a) == len(b)), return alignment result (* for matching, - for non-matching) and longest alignment strech
def getAlignmentMetrics (seq_a, seq_b):
    alignment_result = ""
    count_matching = 0

    for i in range(0, len(seq_a)):
        if seq_a[i] == seq_b[i]:
            alignment_result += "*"
            count_matching += 1
        else:
            alignment_result += "-"

    longest_strech = 0
    current_strech = 0
    in_strech = False

    for i in range(0, len(alignment_result)):

        if alignment_result[i] == "*":
            if in_strech:
                current_strech += 1
            else:
                in_strech = True
                current_strech = 1
        else:
            if in_strech:
                in_strech = False
                if current_strech > longest_strech:
                    longest_strech = current_strech

    if in_strech:
        in_strech = False
        if current_strech > longest_strech:
            longest_strech = current_strech

    return [seq_a, seq_b, alignment_result, count_matching, longest_strech]

# Shifts sequences a and b to each other, calculates algnemt and overlap and returns best match
def getMicrohomologyBaseCount (seq_a, seq_b):

    best_longest_strech = 0
    best_result = []

    for offset in range(len(seq_a)):
        N_offset = ""
        for i in range(0, offset):
            N_offset += "N"

        al_1_results = getAlignmentMetrics(seq_a + N_offset, N_offset + seq_b)
        if al_1_results[4] > best_longest_strech:
            best_longest_strech = al_1_results[4]
            best_result = al_1_results

        al_2_results = getAlignmentMetrics(N_offset + seq_a, seq_b + N_offset)
        if al_2_results[4] > best_longest_strech:
            best_longest_strech = al_2_results[4]
            best_result = al_2_results

    return best_result

# Calculates distance to next LCR
def getDistanceToNextLCR (chrom, pos, LCR_bed_file_path):
    min_distance = 100000000
    min_repeats_name = ""
    min_lcr_dir = "."
    for line in open(LCR_bed_file_path):
        lcr_chrom = line.split("\t")[0]
        lcr_start = int(line.split("\t")[1])
        lcr_end = int(line.split("\t")[2])
        lcr_dir = line.split("\t")[5]

        if chrom == lcr_chrom:
            if pos >= lcr_start and pos <= lcr_end:
                return 0, line.split("\t")[3], lcr_dir
            else:
                dist = np.amin([np.abs(pos - lcr_start), np.abs(pos - lcr_end)])
                if dist < min_distance:
                    min_distance = dist
                    min_repeats_name = line.split("\t")[3]
                    min_lcr_dir = lcr_dir
    return min_distance, min_repeats_name, min_lcr_dir