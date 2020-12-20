#!/bin/sh

# pre1
#bowtie2-inspect /opt/conda/bin/metaphlan_databases/mpa_v20_m200 > /tmp/db_markers/all_markers.fasta

# pre2
# while read line;
# do
# 	extract_markers.py --mpa_pkl /opt/conda/bin/metaphlan_databases/mpa_v20_m200.pkl --ifn_markers /tmp/db_markers/all_markers.fasta --clade $line --ofn_markers /tmp/db_markers/$line\.fasta
# done < /tmp/clades.txt

#extract_markers.py --mpa_pkl /opt/conda/bin/metaphlan_databases/mpa_v20_m200.pkl --ifn_markers /tmp/db_markers/all_markers.fasta --clade g__Aneurinibacillus --ofn_markers /tmp/db_markers/g__Aneurinibacillus.fasta

# pre4
# 1

strain="g__Aneurinibacillus"
wgs="GCF_006539965.1_ASM653996v1_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Anaerotruncus_sp_G3_2012"
wgs="GCF_000403395.2_Anae_bact_G3_V1_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Bacillus_cereus_thuringiensis"
wgs="GCF_000008505.1_ASM850v1_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Clostridium_sporogenes"
wgs="GCF_000960175.1_ASM96017v1_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Desulfotomaculum_ruminis"
wgs="GCF_000215085.1_ASM21508v1_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Enterobacteria_phage_lambda"
wgs="GCA_002745415.1_ASM274541v1_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Enterobacteria_phage_phiX174_sensu_lato"
wgs="GCF_000914915.1_ViralProj240593_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Enterococcus_faecalis"
wgs="GCF_902161805.1_25426_7_320_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Human_endogenous_retrovirus_K"
wgs="GCF_000913595.1_ViralProj222261_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Lactobacillus_johnsonii"
wgs="GCF_003316915.1_ASM331691v1_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Oscillibacter_sp_1_3"
wgs="GCF_000403435.2_Osci_bact_1-3_V1_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

strain="s__Saccharomyces_cerevisiae_killer_virus_M1"
wgs="GCF_000848385.1_ViralProj14678_genomic.fna.gz"

strainphlan.py --ifn_samples /tmp/*.markers --ifn_markers /tmp/db_markers/$strain\.fasta \
--ifn_ref_genomes /tmp/db_markers/$wgs \
 --output_dir /tmp/tree/ --clades $strain --nprocs_main 36 --marker_in_clade 0.5

