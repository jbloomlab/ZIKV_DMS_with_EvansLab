# Summary of steps to recreate PyMOL figures

Run the python script `write_pymol_scripts.py`. This will create three PyMOL scripts that describe the color mapping for structure `5IRE.pdb`. These scripts will map mean site diffsel for escape from the antibodies ZKA64 and ZKA185 (white to red), and will map mutational tolerance as entropy (white to green).
- Open `5IRE.pdb` in PyMOL
- Drag the script into the PyMOL window to color the surface of the dimer and save .png images of top and side views of the structure.
