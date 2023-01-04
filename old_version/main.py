import Filter

if __name__ == '__main__':

    vcf_file_path = "/home/ubuntu/seq/Miltenyi_ONT/sniffles2_final/Miltenyi.SAC.minimap.Y.sniffles.sort.vcf"

    germline_logfile = open("germline.final.logfile.txt", "w")

    logfile = open("logfile.final.txt", "w")

    #

    Filter.reformat_vcf(vcf_file_path, vcf_file_path.replace(".vcf", ".reformat.vcf"))

    Filter.first_stage_filter(vcf_file_path.replace(".vcf", ".reformat.vcf"), vcf_file_path.replace(".vcf", ".reformat.1.filtered.vcf"))

    Filter.annotate_vcf(vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.vcf"), vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.vcf"), "/home/ubuntu/seq/Miltenyi_ONT/mappings/SAC_minimap_Y/Miltenyi.SAC.minimap.Y.bam", logfile)

    Filter.second_stage_filter(vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.vcf"), vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.2.filtered.vcf"))

    Filter.combine_duplicates(vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.2.filtered.pass.vcf"), vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.2.filtered.pass.copies_removed.vcf"))
    
    Filter.germline_filter(vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.2.filtered.pass.copies_removed.vcf"), vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.2.filtered.pass.copies_removed.germline.vcf"), "/home/ubuntu/seq/Miltenyi_ONT/mappings/SAC_minimap_Y/Miltenyi.SAC.minimap.Y.bam", germline_logfile)
    
    Filter.third_stage_filter(vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.2.filtered.pass.copies_removed.germline.vcf"), vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.2.filtered.pass.copies_removed.germline.filtered.vcf"))

    Filter.exp_filter(vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.2.filtered.pass.copies_removed.germline.filtered.pass.vcf"), vcf_file_path.replace(".vcf", ".reformat.1.filtered.pass.annotated.2.filtered.pass.copies_removed.germline.filtered.experimentell.vcf"))
