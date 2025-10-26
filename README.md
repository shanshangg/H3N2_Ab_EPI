# H3N2_Ab_EPI
The data and scripts employed in the research on rapid identification of influenza A virus HA-targeting antibody epitopes.

## Data
1. Hemagglutinin_related_sabdab_summary.tsv
The summary file provided by SAbDab contains multiple properties annotated for each structure. From this file, only entries where the antigen_species was annotated as Influenza A Virus and the antigen_name was annotated as Hemagglutinin were retained. Subsequently, packing complexes were eliminated and antibody redundancies removed using the codes provided by [Ab_interface_mapping](https://github.com/andreasvisbech/Ab_interface_mapping).
2. Sequences_of_antibodies_against_HA.csv
The table provides the amino acid sequences of antibodies, which serve as inputs for antibody structure modeling using ABodyBuilder2.

## Scripts
1. Parses_PDB.py
The script extracts a unique set of antigen and antibody chains from PDB files and saves them as new PDB files.
2. Get_Binding_Site_List.py
The script is used to obtain a list of amino acid residues at the antigen-antibody binding sites.
3. Get_Epitope_Structure.py
The script is used for the batch extraction of epitope structures from antigen-antibody complex structures. It is adapted  from [CE_BLAST/extract_epitope_from_complex.py](https://github.com/baddtongji/CE_BLAST), with dependencies obtained from github repository [CE_BLAST](https://github.com/baddtongji/CE_BLAST).
4. Get_TMscore_among_Epitopes.py
The script is used to batch calculate the structural similarity metric (TM-score) between epitope structures.
5. Get_Chain_Sequences_from_PDB.py
The script is used to extract the amino acid sequences of antibodies and antigens from PDB files.
6. Run_ABodyBuilder2_for_ALL_Antibodies.py
The script is used to perform batch antibody structure modeling using [ABodyBuilder2](https://github.com/oxpig/ImmuneBuilder).