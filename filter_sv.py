import argparse
import prefilter
from time import time, ctime, time_ns
import germline_filter
import os

parser = argparse.ArgumentParser(description='Filters false positive structural variants from a vcf file')
parser.add_argument("--vcf", required=True, help="Path to SV vcf file")
parser.add_argument("--bam", required=True, help="Path to corresponding BAM file")
parser.add_argument("--ref", required=True, help="Path to reference FASTA file")
parser.add_argument("--folder", required=True, help="Path to data folder")

args = parser.parse_args()

data_folder = args.folder
if not data_folder.endswith("/"):
    data_folder += "/"

print("Starting Filter SV at " + ctime(time()))
print("Starting the prefilter module at " + ctime(time()))
time_ns_prefilter_start = time_ns()
total_variants, short_variants, filtered_variants, passing_variants = prefilter.prefilter(args.vcf, args.bam, data_folder)
time_ns_prefilter_end = time_ns()
print("Completing the prefilter module at " + ctime(time()) + ", elapsed time " + str(round((time_ns_prefilter_end - time_ns_prefilter_start) / 1000000000, 3)) + " s")
print("\tTotal variants\t\t" + str(total_variants))
print("\tShort variants\t\t" + str(short_variants))
print("\tFiltered variants\t" + str(filtered_variants))
print("\tPassing variants\t" + str(passing_variants))

print("Starting the germline insertion module at " + ctime(time()))
time_ns_germlinefilter_start = time_ns()
germline_filter.germline_filter (
    args.ref,
    args.vcf.split(".vcf")[0].split("/")[-1] + ".prefilter.vcf",
    args.bam,
    args.folder)
time_ns_germlinefilter_end = time_ns()
print("Completing the germline insertion module at " + ctime(time()) + ", elapsed time " + str(round((time_ns_germlinefilter_end - time_ns_germlinefilter_start) / 1000000000, 3)))

# map to reference
print("Mapping reads and assemblies to the reference genome at " + ctime(time()))
time_ns_mapping_start = time_ns()

os.system("cat " + args.folder + "sv_fastq/*.fastq > tmp.fastq")
os.system("minimap2 -ax map-ont -t 26 " + args.ref + " tmp.fastq -Y | samtools view -b - | samtools sort > " + args.folder + "sv_mappings/sv_reads_to_reference.bam")
os.system("samtools index " + args.folder + "sv_mappings/sv_reads_to_reference.bam")

os.system("cat " + args.folder + "sv_assemblies/*.fasta > tmp.fasta")
os.system("minimap2 -ax map-ont -t 26 " + args.ref + " tmp.fasta -Y | samtools view -b - | samtools sort > " + args.folder + "sv_mappings/sv_assemblies_to_reference.bam")
os.system("samtools index " + args.folder + "sv_mappings/sv_assemblies_to_reference.bam")

os.system("cat " + args.folder + "germline_fastq/*.fastq > tmp.fastq")
os.system("minimap2 -ax map-ont -t 26 " + args.ref + " tmp.fastq -Y | samtools view -b - | samtools sort > " + args.folder + "germline_mappings/germline_reads_to_reference.bam")
os.system("samtools index " + args.folder + "sv_mappings/sv_reads_to_reference.bam")

os.system("cat " + args.folder + "germline_assemblies/*.fasta > tmp.fasta")
os.system("minimap2 -ax map-ont -t 26 " + args.ref + " tmp.fasta -Y | samtools view -b - | samtools sort > " + args.folder + "germline_mappings/germline_assemblies_to_reference.bam")
os.system("samtools index " + args.folder + "germline_mappings/germline_assemblies_to_reference.bam")

os.system("rm tmp.fastq")
os.system("rm tmp.fasta")

time_ns_mapping_end = time_ns()
print("Completing the mapping step at at " + ctime(time()) + ", elapsed time " + str(round((time_ns_mapping_end - time_ns_mapping_start) / 1000000000, 3)))


