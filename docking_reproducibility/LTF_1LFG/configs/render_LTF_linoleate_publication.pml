reinitialize
load ./08_PLIP_PyMOL/LTF_1LFG_Nlobe_linoleate_run02_mode1_complex.pdb, N_complex

hide everything
show cartoon, N_complex and chain A
color gray80, N_complex and chain A
set cartoon_transparency, 0.45

select N_ligand, N_complex and chain L and resn LIA
show sticks, N_ligand
color orange, N_ligand

select N_pocket, byres (N_complex and chain A within 4.0 of N_ligand)
show sticks, N_pocket
color cyan, N_pocket

select N_HIS253, N_complex and chain A and resi 253
show sticks, N_HIS253
color marine, N_HIS253

select N_iron, N_complex and chain A and resn FE
show spheres, N_iron
color magenta, N_iron
set sphere_scale, 0.35, N_iron

distance N_HIS253_contact, (N_ligand and elem O), (N_HIS253 and elem N), 3.5, 2
color red, N_HIS253_contact

label (N_HIS253 and name CA), "HIS253"
label N_iron, "Fe"

set dash_width, 2.5
set dash_gap, 0.3
set dash_length, 0.25
set label_distance_digits, 2
set stick_radius, 0.18
set label_size, 18
set label_color, black
set ray_shadows, 0
set antialias, 2
bg_color white

orient N_ligand
zoom N_ligand, 11

save ./08_PLIP_PyMOL/LTF_1LFG_Nlobe_linoleate_QC.pse
png ./08_PLIP_PyMOL/LTF_1LFG_Nlobe_linoleate_QC.png, 1800, 1400, 300, 1

reinitialize
load ./08_PLIP_PyMOL/LTF_1LFG_Clobe_linoleate_run03_mode1_complex.pdb, C_complex

hide everything
show cartoon, C_complex and chain A
color gray80, C_complex and chain A
set cartoon_transparency, 0.45

select C_ligand, C_complex and chain L and resn LIA
show sticks, C_ligand
color orange, C_ligand

select C_pocket, byres (C_complex and chain A within 4.0 of C_ligand)
show sticks, C_pocket
color cyan, C_pocket

select C_key_residues, C_complex and chain A and resi 529+640
show sticks, C_key_residues
color marine, C_key_residues

select C_iron, C_complex and chain A and resn FE
show spheres, C_iron
color magenta, C_iron
set sphere_scale, 0.35, C_iron

distance C_polar_contacts, (C_ligand and elem O), (C_key_residues and elem N+O), 3.5, 2
color red, C_polar_contacts

label (C_key_residues and name CA), "%s%s" % (resn,resi)
label C_iron, "Fe"

set dash_width, 2.5
set dash_gap, 0.3
set dash_length, 0.25
set label_distance_digits, 2
set stick_radius, 0.18
set label_size, 18
set label_color, black
set ray_shadows, 0
set antialias, 2
bg_color white

orient C_ligand
zoom C_ligand, 11

save ./08_PLIP_PyMOL/LTF_1LFG_Clobe_linoleate_QC.pse
png ./08_PLIP_PyMOL/LTF_1LFG_Clobe_linoleate_QC.png, 1800, 1400, 300, 1

quit
