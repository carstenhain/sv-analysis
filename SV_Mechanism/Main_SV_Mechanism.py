import BuildSVAndTestData

# Padding 25 for RAG heptamer analysis
"""
BuildSVAndTestData.build_BED_around_SV_breakpoints("C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WGS\\SV\\CD3Iso.Manuell.Truth_Set.plus.Add_In.templates.manuell_kontrolliert.Final.Sortiert.bed", "..\\Resources\\SV.breakpoints.padding.25.bed", True, 25)
BuildSVAndTestData.build_BED_around_random_breakpoints(200, "..\\Resources\\Random.breakpoints.padding.25.bed", 25)
"""
# Padding 12 for microhomology analysis
"""
BuildSVAndTestData.build_BED_around_SV_breakpoints("C:\\Users\\Carst\\Documents\\Promotion\\MF_Verlaufspatient\\WGS\\SV\\CD3Iso.Manuell.Truth_Set.plus.Add_In.templates.manuell_kontrolliert.Final.Sortiert.bed", "..\\Resources\\SV.breakpoints.padding.12.bed", True, 12)
BuildSVAndTestData.build_BED_around_random_breakpoints(200, "..\\Resources\\Random.breakpoints.padding.12.bed", 12)
"""
# bedtools getfasta -fi hg38.fasta -bed padded.bed > padded.fasta
