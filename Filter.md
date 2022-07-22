# Filtering of structural variants (SVs)
All cutoffs are adjusted for this dataset
## Coverage filter I
- Filtering SVs in regions with abnormal high coverage
  - Sniffles2 annotation **DV > 50**
  - Often Sniffles annotation are not accurate or the SV breakpoint is 1-2 bp next to a region with abnormal high coverage. Therefore, mean coverage in a 20 bp segment around each breakpoint is calculate using *samtools bedcov*. SVs with breakpoints in a high coverage segment (**Coverage > 100**) are filtered out.
- Filtering reads with less than 2 or more than 15 supporting reads (**DV < 2 or DV > 15**)
## SV length filter
- Filtering SVs with a length smaller than 5 kb
## PCR duplicate filter
- Assembly of SV supporting reads with lamassemble [(Frith et al., 2021)](https://pubmed.ncbi.nlm.nih.gov/33289891/)
- Mapping of SV supporting reads on the assembly using minimap2
- 