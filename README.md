# sv-analysis
Scripts for SV analysis - including filtering, refining, plotting
## Files
### General
| File Name | File Type | Description |
| ----------- | ----------- | ----------- |
| ChromosomeMetrics.py | Python | Methods for GRCh38 chromosomes, length, positions, ... |
| LoadData.py | Python | Methods for loading different data types from file |
### Resources
| File Name | File Type | Description |
| ----------- | ----------- | ----------- |
| karyotype.human.hg38.txt | Text | Karyotype file (staining) for hg38 from circos 0.69.9, chromosomes as hs1 instead of chr1 |
| RAG.heptamer.meme | Text | Motif file for MEME build from TCR RAG heptamers |
### SV_Mechanism
| File Name | File Type | Description |
| ----------- | ----------- | ----------- |
| BuildSVAndTestData.py | Python | Methods for building files for, and analysing potential SV mechanism. Including RAG enrichment (requires manual FIMO step), LCR enrichment and microhomology  |
| Plot.py | Python | Plot gathered data |
