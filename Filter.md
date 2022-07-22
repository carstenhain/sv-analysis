# Filtering of structural variants (SVs)
## Coverage filter I
- Filtering SVs in regions with abnormal high coverage (DV > 50)
- Filtering reads with less than 2 or more than 15 supporting reads (adjusted for this dataset)
## SV length filter
- Filtering SVs with a length smaller than 5 kb
## PCR duplicate filter
- Assembly of SV supporting reads with lamassemble [(Frith et al., 2021)](https://pubmed.ncbi.nlm.nih.gov/33289891/)
