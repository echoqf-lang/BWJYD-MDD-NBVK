reinitialize
load ./LTF_1LFG/01_receptor_clean/LTF_1LFG_chainA_protein_repaired.pdb, protein
load ./LTF_1LFG/00_raw_inputs/1LFG_raw.pdb, raw_1LFG
select iron_components, raw_1LFG and chain A and (resn FE or resn CO3) and resi 701+702+703+704
create iron_and_carbonate, iron_components
create LTF_1LFG_diferric_receptor, protein or iron_and_carbonate
sort
save ./LTF_1LFG/01_receptor_clean/LTF_1LFG_chainA_2Fe_2CO3_repaired.pdb, LTF_1LFG_diferric_receptor
quit
