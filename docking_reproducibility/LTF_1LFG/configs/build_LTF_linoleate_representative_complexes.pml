reinitialize

load ./01_receptor_clean/LTF_1LFG_chainA_2Fe_2CO3_repaired.pdb, receptor

load ./05_Nlobe_docking/linoleate_Nlobe_run02_poses.sdf, Nlobe_poses
create linoleate_N_top, Nlobe_poses, 1, 1
alter linoleate_N_top, resn="LIA"
alter linoleate_N_top, chain="L"
alter linoleate_N_top, resi="1"
alter linoleate_N_top, segi=""
sort
create LTF_Nlobe_linoleate_complex, receptor or linoleate_N_top
save ./08_PLIP_PyMOL/LTF_1LFG_Nlobe_linoleate_run02_mode1_complex.pdb, LTF_Nlobe_linoleate_complex

load ./06_Clobe_docking/linoleate_Clobe_run03_poses.sdf, Clobe_poses
create linoleate_C_top, Clobe_poses, 1, 1
alter linoleate_C_top, resn="LIA"
alter linoleate_C_top, chain="L"
alter linoleate_C_top, resi="1"
alter linoleate_C_top, segi=""
sort
create LTF_Clobe_linoleate_complex, receptor or linoleate_C_top
save ./08_PLIP_PyMOL/LTF_1LFG_Clobe_linoleate_run03_mode1_complex.pdb, LTF_Clobe_linoleate_complex

quit
