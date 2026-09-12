reinitialize
load ./LCN2_1QQS/01_receptor_clean/LCN2_1QQS_chainA_sidechains_repaired.pdb, receptor
load ./LCN2_1QQS/08_linoleic_docking/linoleate_run03_poses.sdf, linoleate_poses
create linoleate_top, linoleate_poses, 1, 1
alter linoleate_top, resn="LIA"
alter linoleate_top, chain="L"
alter linoleate_top, resi="1"
alter linoleate_top, segi=""
sort
create LCN2_linoleate_complex, receptor or linoleate_top
save ./LCN2_1QQS/09_PLIP_PyMOL/LCN2_1QQS_linoleate_run03_mode1_complex.pdb, LCN2_linoleate_complex
quit
