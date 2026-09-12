reinitialize
load ./08_PLIP_PyMOL/LTF_1LFG_Nlobe_linoleate_run02_mode1_complex.pdb, N_complex

hide everything
show cartoon, N_complex and chain A
color gray80, N_complex and chain A
set cartoon_transparency, 0.55

select N_ligand, N_complex and chain L and resn LIA
show sticks, N_ligand
color orange, N_ligand

select N_hydrophobic, N_complex and chain A and resi 11+60+63+183+190+297
show sticks, N_hydrophobic
color cyan, N_hydrophobic

select N_HIS253, N_complex and chain A and resi 253
show sticks, N_HIS253
color marine, N_HIS253

select N_FE701, N_complex and chain A and resn FE and resi 701
show spheres, N_FE701
color magenta, N_FE701
set sphere_scale, 0.35, N_FE701

distance N_contact, (N_ligand and elem O), (N_HIS253 and elem N), 3.5, 2
color red, N_contact
hide labels, N_contact

label (N_HIS253 and name CA), "HIS253"
set label_position, [2.0, 1.0, 0.0], (N_HIS253 and name CA)
label N_FE701, "Fe701"
set label_position, [-1.5, 1.5, 0.0], N_FE701

set dash_width, 2.5
set dash_gap, 0.3
set dash_length, 0.25
set stick_radius, 0.18
set label_size, 16
set label_color, black
set ray_shadows, 0
set antialias, 2
bg_color white

orient N_ligand
zoom N_ligand, 9

save ./08_PLIP_PyMOL/LTF_1LFG_Nlobe_linoleate_publication.pse
png ./08_PLIP_PyMOL/LTF_1LFG_Nlobe_linoleate_publication.png, 1800, 1400, 300, 1

reinitialize
load ./08_PLIP_PyMOL/LTF_1LFG_Clobe_linoleate_run03_mode1_complex.pdb, C_complex

hide everything
show cartoon, C_complex and chain A
color gray80, C_complex and chain A
set cartoon_transparency, 0.55

select C_hydrophobic, C_complex and chain A and resi 350+392+398+526
show sticks, C_hydrophobic
color cyan, C_hydrophobic

select C_key_residues, C_complex and chain A and resi 529+640
show sticks, C_key_residues
color marine, C_key_residues

select C_ligand, C_complex and chain L and resn LIA
show sticks, C_ligand
color orange, C_ligand

select C_FE702, C_complex and chain A and resn FE and resi 702
show spheres, C_FE702
color magenta, C_FE702
set sphere_scale, 0.35, C_FE702

distance C_contacts, (C_ligand and elem O), (C_key_residues and elem N+O), 3.5, 2
color red, C_contacts
hide labels, C_contacts

label (C_key_residues and resi 529 and name CA), "THR529"
set label_position, [2.5, 1.5, 0.0], (C_key_residues and resi 529 and name CA)

label (C_key_residues and resi 640 and name CA), "ASN640"
set label_position, [-2.5, -1.5, 0.0], (C_key_residues and resi 640 and name CA)

label C_FE702, "Fe702"
set label_position, [-1.5, 1.5, 0.0], C_FE702

set dash_width, 2.5
set dash_gap, 0.3
set dash_length, 0.25
set stick_radius, 0.18
set label_size, 16
set label_color, black
set ray_shadows, 0
set antialias, 2
bg_color white

orient C_ligand
zoom C_ligand, 9

save ./08_PLIP_PyMOL/LTF_1LFG_Clobe_linoleate_publication.pse
png ./08_PLIP_PyMOL/LTF_1LFG_Clobe_linoleate_publication.png, 1800, 1400, 300, 1

quit
