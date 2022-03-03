def load_sv_bed (sv_bed_file_path):

    sv_dict = {}

    for line in open(sv_bed_file_path):
        chrom = line.rstrip().split("\t")[0]
        start = line.rstrip().split("\t")[1]
        end = line.rstrip().split("\t")[2]
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
    for sv in sv_dict:
        print(sv)
        for item in sv_dict[sv]:
            print(item)
        print("-----------------")
    """

