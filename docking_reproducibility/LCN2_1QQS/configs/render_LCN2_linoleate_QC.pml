reinitialize
load ./LCN2_1QQS/09_PLIP_PyMOL/LCN2_1QQS_linoleate_run03_mode1_complex.pdb, complex

hide everything
show cartoon, chain A
color gray80, chain A
set cartoon_transparency, 0.45

select linoleate, chain L and resn LIA
show sticks, linoleate
color orange, linoleate

select pocket, byres (chain A within 4.0 of linoleate)
show sticks, pocket
color cyan, pocket

select key_residues, chain A and resi 81+134
show sticks, key_residues
color marine, key_residues

distance polar_contacts, (linoleate and elem O), (key_residues and elem N), 3.5, 2
color red, polar_contacts
set dash_width, 2.5
set dash_gap, 0.3
set dash_length, 0.25
set label_distance_digits, 2

label (chain A and resi 81+134 and name CA), "%s%s" % (resn,resi)

set stick_radius, 0.18
set label_size, 18
set label_color, black
set ray_shadows, 0
set antialias, 2
bg_color white

orient linoleate
zoom linoleate, 10

save ./LCN2_1QQS/09_PLIP_PyMOL/LCN2_1QQS_linoleate_QC.pse
png ./LCN2_1QQS/09_PLIP_PyMOL/LCN2_1QQS_linoleate_QC.png, 1800, 1400, 300, 1
quit
