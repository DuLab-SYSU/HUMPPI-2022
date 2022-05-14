Biomart_geneid.txt # ID mapping from biomart Ensembl genome browser 100 (2020-05-27)
Interaction_detection_method.txt # The PPI interaction detection methods (binary or non-binary) https://www.ebi.ac.uk/ols/ontologies/MI/ (2020-05-31)
protein_coding_genes.txt # 17,402 protein-coding genes from human ORFeome v9.1 (Luck K, et al. 2020) (2020-06-09)
human_name_age.txt # human gene name and phylogenetic branch, gene have two age and the older is retained (2020-08-17)
ID2uniprot.txt # ID mapping from RefSeq protein/EMBL/GenBank/DDBJ/GeneID/UniProtKB to UniProtKB, https://www.uniprot.org (2020-06-26)
essential_genes.txt # essential gene list from Blomen, et al. (2015) and Wang, et al. (2015) (2021-01-04)
virus_description.txt # the viruses information from NCBI, UniProt and Virus-Host DB (2021-09-18)
GTEx_Analysis_gene_tissue-specific_expression.gct # tissue-specific expression for protein-coding genes (2021-03-05)

cd PPIs_from_databases/
BIOGRID-ORGANISM-Homo_sapiens-3.5.182.mitab.txt.gz # https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.5.182/BIOGRID-ALL-3.5.182.mitab.zip (2020-03-20) 
matrixdb_human.tab.gz # http://matrixdb.univ-lyon1.fr (2020-04-22)
Hsapi20170205.txt.gz # https://dip.doe-mbi.ucla.edu/dip/Download.cgi?SM=7&TX=9606 (2020-03-17)                               
MINT_human.gz # http://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/current/search/query/species:human (2020-03-17)
intact_human.txt.gz # ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip (2020-03-18)

cd ../IRGs/
all_IRGs.txt # the unique gene list from ImmPort, InnateDB, Immunome database and PathCards (2020-07-18)
Immune-related_process.txt # unique same pathway from different databases (Supplementary table 2) (2021-02-04)

cd ../MSigDB7.4
c5.bp.v7.4.symbols_GO.ID.txt.gz # The gene ontology terms were obtained from the c5 category of Molecular Signature Database (MSigDB v7.4) (2021-08-10) 

cd ../VHI_from_databases
BioGRID_human.txt.gz, intact_virus.txt.gz, intact_human-others.txt.gz # see information above
hpidb2.mitab.txt.gz # HPIDB 3.0, https://hpidb.igbb.msstate.edu/about.html#cites, export LC_COLLATE='C' export LC_CTYPE='C' cat hpidb2.mitab.txt | cut -f1,2,7,9,10,11,12, solve the problem: 'stdin: Illegal byte sequence' (2020-06-18)
virhostnet.txt.gz # VirHostNet 2.0, http://virhostnet.prabi.fr (2020-06-16) 
From_literatures/ # Coronavirus-human PPI from 3 literatures (2021-09-17)


