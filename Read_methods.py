import pysam
import numpy as np
import os
import subprocess
import pysam
import matplotlib.pyplot as plt

def collect_reads (read_names, chrom, pos, search_area, alignment_file, logfile_fo):
    logfile_fo.write("\t\t\tLooking at interval " + chrom + ":" + str(np.amax([0, pos - search_area])) + "-" + str(pos + search_area) + "\n")
    supporting_alignments = []
    for read in alignment_file.fetch(chrom, np.amax([0, pos - search_area]), pos + search_area):
        if read.query_name in read_names:
            read_already_saved = False
            for r in supporting_alignments:
                if r.query_name == read.query_name:
                    read_already_saved = True
            if not read_already_saved:
                supporting_alignments.append(read)
    if len(supporting_alignments) == len(read_names):
        logfile_fo.write("\t\t\tFound all supporting reads\n")
    else:
        logfile_fo.write("\t\t\t" + str(len(read_names) - len(supporting_alignments)) + " reads missing\n")
        reads = read_names
        for s in supporting_alignments:
            reads.remove(s.query_name)
        logfile_fo.write("\n".join(["\t\t\t\t" + str(s) for s in reads]) + "\n")
    return supporting_alignments

def write_fastq (read, fastqfile_fo):
    fastqfile_fo.write("@%s\n" % read.qname)
    fastqfile_fo.write("%s\n" % read.seq)
    fastqfile_fo.write("+\n")
    fastqfile_fo.write("%s\n" % read.qual)
    fastqfile_fo.flush()

def write_fastq_original (read, fastqfile_fo):
    if not read.is_reverse:
        write_fastq(read, fastqfile_fo)
    else:
        fastqfile_fo.write("@%s\n" % read.qname)
        fastqfile_fo.write("%s\n" % read.get_forward_sequence())
        fastqfile_fo.write("+\n")
        for_qualities = []
        rev_qualities = read.qual
        for i in range(0, len(rev_qualities)):
            for_qualities.append(rev_qualities[len(rev_qualities) - i - 1])
        fastqfile_fo.write("%s\n" % "".join(for_qualities))
        fastqfile_fo.flush()

def get_reads_mean_MQ (reads):
    mq = []
    for read in reads:
        mq.append(read.mapping_quality)
    return np.mean(mq)

def get_read_mean_mismatch_fraction (reads):
    mismatch_fractions = []
    for read in reads:
        num_mismatches = 0
        for tup in read.get_aligned_pairs(with_seq=True):
            if tup[2] in ["a", "c", "g", "t"]:
                num_mismatches += 1
        mismatch_fractions.append(num_mismatches / read.reference_length)
    return np.mean(mismatch_fractions)

def read_single_fasta (fasta_file_path):
    fasta_name = ""
    fasta_sequence = ""
    for line in open(fasta_file_path):
        if line.startswith(">"):
            if fasta_name == "":
                fasta_name = line[1:].rstrip()
            else:
                break
        else:
            fasta_sequence += line.rstrip()
    return fasta_name, fasta_sequence

def assemble_reads (path_to_fastq, path_to_lamassemble, path_to_lamassemble_train_file, assembly_name, savepath_for_assembly_fasta, logfile_fo):
    # assemble reads
    os.system(path_to_lamassemble + " -P 8 --all -n " + assembly_name + " " + path_to_lamassemble_train_file + " " + path_to_fastq + " > " + savepath_for_assembly_fasta)

    logfile_fo.write("\t\tStarting assembly\n")

    # QC if assembly produced contig
    a_name, a_seq = read_single_fasta(savepath_for_assembly_fasta)
    logfile_fo.write("\t\t\tAssembly produced a contig with length " + str(len(a_seq)) + "\n")

def detect_duplicates (sv_record, supporting_reads, path_to_fastq, path_to_assembly_fasta, path_to_minimap, logfile_fo):

    logfile_fo.write("\t\tMapping reads into mappings/" + sv_record.ID + ".bam\n")
    print()
    os.system(path_to_minimap + " -ax map-ont -Y -t 10 " + path_to_assembly_fasta + " " + path_to_fastq + " | samtools view -b - | samtools sort - > mappings/" + sv_record.ID + ".bam")
    os.system("samtools index mappings/" + sv_record.ID + ".bam")

    samfile = pysam.AlignmentFile("mappings/" + sv_record.ID + ".bam", "rb")
    reads_on_assembly = []
    for read in samfile.fetch():
        reads_on_assembly.append(read)

    plt.subplots()
    a_name, a_seq = read_single_fasta(path_to_assembly_fasta)
    plt.xlim(0, len(a_seq))
    plt.xlabel("Position")

    i = 0
    for read in samfile.fetch():
        plt.plot([read.reference_start, read.reference_end], [i, i])
        i += 1

    plt.ylim(-1, i)
    plt.savefig("mapping_plots/" + sv_record.ID + ".png")
    plt.close()

    read_start_distances = np.zeros((len(reads_on_assembly), len(reads_on_assembly)))
    read_end_distances = np.zeros((len(reads_on_assembly), len(reads_on_assembly)))
    is_duplicate = np.zeros((len(reads_on_assembly), len(reads_on_assembly)))

    for i in range(0, len(reads_on_assembly)):
        for j in range(0, len(reads_on_assembly)):
            if j != i:
                logfile_fo.write("\t\t\tComparing reads " + reads_on_assembly[i].query_name + " with " + reads_on_assembly[j].query_name + "\n")
                logfile_fo.flush()

                # Calculate reads start and end distances

                distance_start = np.abs(reads_on_assembly[i].reference_start - reads_on_assembly[j].reference_start)
                distance_end = np.abs(reads_on_assembly[i].reference_end - reads_on_assembly[j].reference_end)

                # check for same library

                lib_i = ""
                lib_j = ""

                for read in supporting_reads:
                    if read.query_name == reads_on_assembly[i].query_name:
                        lib_i = read.get_tag("RG")
                    if read.query_name == reads_on_assembly[j].query_name:
                        lib_j = read.get_tag("RG")

                if lib_i != lib_j:
                    distance_start = 1000000
                    distance_end = 1000000
                    logfile_fo.write("\t\t\t\tReads come from different libraries\n")
                else:
                    logfile_fo.write("\t\t\t\tDistance start: " + str(distance_start) + "\n")
                    logfile_fo.write("\t\t\t\tDistance end: " + str(distance_end) + "\n")

                #Calculate if reads are duplicates based on start and end distance

                read_start_distances[i, j] = distance_start
                read_end_distances[i, j] = distance_end
                if distance_start + distance_end < 150:
                    is_duplicate[i, j] = 1000
    """
    logfile_fo.write("Distance start\n")
    logfile_fo.write(str(read_start_distances) + "\n")
    logfile_fo.write("Distance end\n")
    logfile_fo.write(str(read_end_distances) + "\n")
    logfile_fo.write("Is duplicate\n")
    logfile_fo.write(str(is_duplicate) + "\n")
    """

    plt.subplot(131)
    max = np.amax([1000, np.amax(read_start_distances), np.amax(read_end_distances)])
    plt.imshow(read_start_distances, vmin=0, vmax=200)
    plt.xlabel("Distance start")
    plt.subplot(132)
    plt.imshow(read_end_distances, vmin=0, vmax=200)
    plt.xlabel("Distance end")
    plt.subplot(133)
    plt.imshow(is_duplicate, vmin=0, vmax=200)
    plt.xlabel("Is duplicate")
    plt.colorbar()
    plt.savefig("duplicate_plots/" + sv_record.ID + ".png")
    plt.close()

    rmdup_rnames = []

    # list indicating non duplicate state of read: 0 no duplicate or best read out of a duplicate group, 1: duplicate
    read_duplicate_list = np.zeros(len(reads_on_assembly))
    # list with mean base call qualities of reads
    read_bc_quality = np.zeros(len(reads_on_assembly))
    for i in range(0, len(reads_on_assembly)):
        read_bc_quality[i] = np.mean(reads_on_assembly[i].query_qualities)

    # cycle all read pairs if they are not already listed as duplicates
    # look if there are other reads indicated as duplicates and flag them as duplicates
    for i in range(0, len(reads_on_assembly)):
        if read_duplicate_list[i] == 0:
            duplicates = [i]
            for j in range(0, len(reads_on_assembly)):
                if read_duplicate_list[j] == 0:
                    if is_duplicate[i, j] > 0:
                        duplicates.append(j)
            # duplicates found -> look for read with highest mean base call quality, retain this read and flag the other as duplicates
            if len(duplicates) > 1:
                max_mean_bc_quality = 0
                max_z = -1
                # find read with highest base call quality
                for z in range(0, len(duplicates)):
                    if max_mean_bc_quality < read_bc_quality[duplicates[z]]:
                        max_mean_bc_quality = read_bc_quality[duplicates[z]]
                        max_z = z
                for z in range(0, len(duplicates)):
                    if z != max_z:
                        read_duplicate_list[duplicates[z]] = 1


    logfile_fo.write("\t\t\t" + str(np.sum(read_duplicate_list)) + " duplicates found\n")
    logfile_fo.write("\t\t\tDuplicate list\n")
    logfile_fo.write("\t\t\t" + str(read_duplicate_list) + "\n")

    # add all reads not flagged as duplicates to return output

    rmdup_rnames = []

    for i in range(0, len(reads_on_assembly)):
        if read_duplicate_list[i] == 0:
            rmdup_rnames.append(reads_on_assembly[i].query_name)


    # For SV with 2 read (DV and RMDV) from the same library

    dist_start = 10000
    dist_end = 10000

    if len(rmdup_rnames) == 2 and len(reads_on_assembly) == 2:

        lib_i = ""
        lib_j = ""

        for read in supporting_reads:
            if read.query_name == rmdup_rnames[0]:
                lib_i = read.get_tag("RG")
            if read.query_name == rmdup_rnames[1]:
                lib_j = read.get_tag("RG")

        if lib_i == lib_j:
            dist_start = np.abs(reads_on_assembly[0].reference_start - reads_on_assembly[1].reference_start)
            dist_end = np.abs(reads_on_assembly[0].reference_end - reads_on_assembly[1].reference_end)
            logfile_fo.write("\t\t\t\tReturning data for DIST_START and DIST_END filtering\n")
            logfile_fo.write("\t\t\t\tDistance start: " + str(dist_start) + "\n")
            logfile_fo.write("\t\t\t\tDistance end: " + str(dist_end) + "\n")

    return rmdup_rnames, dist_start, dist_end

def get_end_std_dev (reads, logfile_fo):
    left = []
    right = []

    left_clipping = []
    right_clipping = []

    for read in reads:
        left.append(read.reference_start)
        right.append(read.reference_end)

        if read.cigartuples[0][0] == 4:
            left_clipping.append(read.cigartuples[0][1])
        else:
            left_clipping.append(0)

        if read.cigartuples[-1][0] == 4:
            right_clipping.append(read.cigartuples[-1][1])
        else:
            right_clipping.append(0)

    if np.mean(left_clipping) > np.mean(right_clipping):
        return int(np.mean(left)), np.std(left)
    else:
        return int(np.mean(right)), np.std(right)



