reinitialize
load ./ARG1_3KV2/00_raw_inputs/3KV2_raw.pdb, raw_3KV2
select receptor_selection, (raw_3KV2 and polymer.protein and chain A) or (raw_3KV2 and resn MN and chain A and resi 514+515)
create ARG1_3KV2_chainA_2Mn, receptor_selection
remove ARG1_3KV2_chainA_2Mn and solvent
sort
save ./ARG1_3KV2/01_receptor_clean/ARG1_3KV2_chainA_2Mn_clean.pdb, ARG1_3KV2_chainA_2Mn
quit
