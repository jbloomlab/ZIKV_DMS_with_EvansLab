# Summary of steps to recreate Pymol figures

Run the python script `write_pymol_scripts.py`. This will create four pymol scripts that describe the color mapping for structure `5IRE.pdb`. These scripts will map mean site diffsel for escape from the antibodies ZKA64 and ZKA185, and will map mutational tolerance as entropy and neffective. To create colored dimer images, carry out the following steps:
- Open `5IRE.pdb` in pymol
- To isolate the two chains in the dimer and create a new object called "dimer", type the command `extract dimer, chain C+E`
- To show surface view type `show surface`
- To hide the extraneous parts of the structure type `hide everything, 5IRE`
- To map coloring on the structure, drag the python script with the color assignments into the Pymol window.
- To show a top-view of the structure, copy and paste the following into Pymol:
```
set_view (\
    -0.352205604,    0.911581576,   -0.212054998,\
     0.782968819,    0.162855446,   -0.600361347,\
    -0.512743294,   -0.377484351,   -0.771101713,\
     0.000000000,    0.000000000, -447.172607422,\
  -115.573463440, -110.383399963, -125.525375366,\
   352.554290771,  541.790893555,  -20.000000000 )
```
- To show a side-view of the structure, copy and paste the following into Pymol:
```
set_view (\
    -0.375719965,   -0.311803222,   -0.872701287,\
     0.787538230,   -0.603787839,   -0.123329692,\
    -0.488471270,   -0.733628869,    0.472416550,\
     0.000000000,    0.000000000, -447.172607422,\
  -115.573463440, -110.383399963, -125.525375366,\
  -21650.402343750, 22544.750000000,  -20.000000000 )
```
- To save an image, click File > Export Image As > PNG... and select 'ray trace with transparent background'
