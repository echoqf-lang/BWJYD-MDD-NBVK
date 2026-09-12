reinitialize
load ./LCN2_1QQS/00_raw_inputs/1QQS_raw.pdb, raw_1QQS
select receptor_selection, raw_1QQS and polymer.protein and chain A
create LCN2_1QQS_chainA, receptor_selection
remove LCN2_1QQS_chainA and solvent
sort
save ./LCN2_1QQS/01_receptor_clean/LCN2_1QQS_chainA_clean.pdb, LCN2_1QQS_chainA
quit

