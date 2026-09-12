reinitialize
load ./ARG1_3KV2/01_receptor_clean/ARG1_3KV2_chainA_2Mn_clean.pdb, receptor
load ./ARG1_3KV2/08_puerarin_docking/puerarin_run01_poses.sdf, puerarin_poses
create puerarin_top, puerarin_poses, 1, 1
alter puerarin_top, resn="PUE"
alter puerarin_top, chain="L"
alter puerarin_top, resi="1"
alter puerarin_top, segi=""
sort
create ARG1_puerarin_complex, receptor or puerarin_top
save ./ARG1_3KV2/09_PLIP_PyMOL/ARG1_3KV2_puerarin_run01_mode1_complex.pdb, ARG1_puerarin_complex
quit
