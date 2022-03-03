import numpy as np
import random

# Returns all canonical human chromosomes as array
def get_chromosomes():
    return ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
            "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

# Returns all canonical human chromosomes minus chrX and chrY as array
def get_chromosomes_woXY():
    return ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
            "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

# Returns all canonical human chromosome with their length as dictionary
def get_chromosome_lengths():
    return {
        "chr1": 248956422,
        "chr2": 242193529,
        "chr3": 198295559,
        "chr4": 190214555,
        "chr5": 181538259,
        "chr6": 170805979,
        "chr7": 159345973,
        "chr8": 145138636,
        "chr9": 138394717,
        "chr10": 133797422,
        "chr11": 135086622,
        "chr12": 133275309,
        "chr13": 114364328,
        "chr14": 107043718,
        "chr15": 101991189,
        "chr16": 90338345,
        "chr17": 83257441,
        "chr18": 80373285,
        "chr19": 58617616,
        "chr20": 64444167,
        "chr21": 46709983,
        "chr22": 50818468,
        "chrX": 156040895,
        "chrY": 57227415
    }

# Returns [start, end] for a chromosomal arm
def get_arm_positions(chrom, arm):
    arm_positions = {
        "chr1_p": [0, 123400000],
        "chr1_q": [123400001, 248956422],

        "chr2_p": [0, 93900000],
        "chr2_q": [93900001, 242193529],

        "chr3_p": [0, 90900000],
        "chr3_q": [90900001, 198295559],

        "chr4_p": [0, 50000000],
        "chr4_q": [50000001, 190214555],

        "chr5_p": [0, 48800000],
        "chr5_q": [48800001, 181538259],

        "chr6_p": [0, 59800000],
        "chr6_q": [59800001, 170805979],

        "chr7_p": [0, 60100000],
        "chr7_q": [60100001, 159345973],

        "chr8_p": [0, 45200000],
        "chr8_q": [45200001, 145138636],

        "chr9_p": [0, 43000000],
        "chr9_q": [43000001, 138394717],

        "chr10_p": [0, 39800000],
        "chr10_q": [39800001, 133797422],

        "chr11_p": [0, 53400000],
        "chr11_q": [53400001, 135086622],

        "chr12_p": [0, 35500000],
        "chr12_q": [35500001, 133275309],

        "chr13_p": [0, 17700000],
        "chr13_q": [17700001, 114364328],

        "chr14_p": [0, 17200000],
        "chr14_q": [17200001, 107043718],

        "chr15_p": [0, 19000000],
        "chr15_q": [19000001, 101991189],

        "chr16_p": [0, 36800000],
        "chr16_q": [36800001, 90338345],

        "chr17_p": [0, 25100000],
        "chr17_q": [27300001, 83200000],

        "chr18_p": [0, 18500000],
        "chr18_q": [18500001, 80373285],

        "chr19_p": [0, 26200000],
        "chr19_q": [26200001, 58617616],

        "chr20_p": [0, 28100000],
        "chr20_q": [28100001, 64444167],

        "chr21_p": [0, 12000000],
        "chr21_q": [12000001, 46709983],

        "chr22_p": [0, 15000000],
        "chr22_q": [15000001, 50818468]
    }

    return arm_positions[chrom + "_" + arm]

# Returns Giemsa stain results for a position
def get_stain (chrom, pos, circos_hg38_karyotype_file_path):
    for line in open(circos_hg38_karyotype_file_path):
        if line.startswith("band"):
            band_chrom = line.split(" ")[1].replace("hs", "chr")
            band_start = int(line.split(" ")[4])
            band_end = int(line.split(" ")[5])
            band_stain = line.rstrip().split(" ")[-1]

            if chrom == band_chrom and pos >= band_start and pos <= band_end:
                return band_stain
    return "FAIL"

#Return complete GRCh38 length
def get_genome_length():
    return 3088269832

# Returns absolute position for a chromosome + position combination as int
def get_absolute_position(chrom, pos):
    abs_pos = np.int64(pos)

    for chromosome in get_chromosomes():

        if chromosome == chrom:
            break

        abs_pos += np.int64(get_chromosome_lengths()[chromosome])

    return abs_pos

# Creates a random position on the the genome outside of centromere and variable regions
def get_random_position(circos_hg38_karyotype_file_path):

    for i in range(0, 1000):
        rand_abs_position = random.randrange(1, get_genome_length())

        chrom = ""
        pos = 0

        test_abs_pos = 0

        for chromosome in get_chromosomes():
            test_abs_pos += np.int64(get_chromosome_lengths()[chromosome])
            if test_abs_pos > rand_abs_position:
                chrom = chromosome
                pos = rand_abs_position - test_abs_pos + get_chromosome_lengths()[chromosome]
                break

        print(rand_abs_position)
        print(chrom)
        print(pos)

        stain = get_stain(chrom, pos, circos_hg38_karyotype_file_path)

        if not stain in ["acen", "gvar", "stalk", "FAIL"]:
            return chrom, pos

    raise GeneratorExit("No position outside of centromeres or variable regions found after 1000 tries")
