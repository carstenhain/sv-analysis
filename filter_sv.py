import argparse
import prefilter
from time import time, ctime, time_ns

parser = argparse.ArgumentParser(description='Filters false positive structural variants from a vcf file')
parser.add_argument("--vcf", required=True, help="Path to SV vcf file")
parser.add_argument("--bam", required=True, help="Path to corresponding BAM file")

args = parser.parse_args()

print("Starting Filter SV at " + ctime(time()))
print("Starting the prefilter module at " + ctime(time()))
time_ns_prefilter_start = time_ns()
total_variants, short_variants, filtered_variants, passing_variants = prefilter.prefilter(args.vcf, args.bam)
time_ns_prefilter_end = time_ns()
print("Completing the prefilter module at " + ctime(time()) + ", elapsed time " + str(round((time_ns_prefilter_end - time_ns_prefilter_start) / 1000000000, 3)))
print("\tTotal variants\t\t" + str(total_variants))
print("\tShort variants\t\t" + str(short_variants))
print("\tFiltered variants\t" + str(filtered_variants))
print("\tPassing variants\t" + str(passing_variants))
