reinitialize
load ./ARG1_3KV2/09_PLIP_PyMOL/ARG1_3KV2_puerarin_run01_mode1_complex.pdb, complex

hide everything
show cartoon, chain A
color gray80, chain A

select puerarin, chain L and resn PUE
show sticks, puerarin
color tv_yellow, puerarin

select manganese, chain A and resn MN
show spheres, manganese
color magenta, manganese
set sphere_scale, 0.35, manganese

select pocket, byres (chain A within 4.0 of puerarin)
show sticks, pocket
color cyan, pocket

select key_residues, chain A and resi 130+137+139
show sticks, key_residues
color marine, key_residues

distance polar_contacts, (puerarin and elem O), (key_residues and elem N+O), 3.5, 2
color red, polar_contacts
set dash_width, 2.5
set dash_gap, 0.3
set dash_length, 0.25
set label_distance_digits, 2

label (chain A and resi 130+137+139 and name CA), "%s%s" % (resn,resi)
label manganese, "Mn"

set cartoon_transparency, 0.35
set stick_radius, 0.18
set label_size, 18
set label_color, black
set ray_shadows, 0
set antialias, 2
bg_color white

orient puerarin
zoom puerarin, 11

save ./ARG1_3KV2/09_PLIP_PyMOL/ARG1_3KV2_puerarin_QC.pse
png ./ARG1_3KV2/09_PLIP_PyMOL/ARG1_3KV2_puerarin_QC.png, 1800, 1400, 300, 1
quit
