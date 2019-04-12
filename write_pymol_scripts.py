'''Script modified by Danny Lawrence 2019 to write pymol scripts that will color structure 5IRE of ZIKV E protein. '''

import pandas as pd 
from colour import Color
import os
import dms_tools2


#lets try to map to structure
diffseldir = './results/diffsel/'
pymoldir = './pymol_scripts/'
if not os.path.isdir(pymoldir):
    os.mkdir(pymoldir)

entropy_file = './results/prefs/unscaled_prefs_entropy.csv'
if not os.path.isfile(entropy_file):
    charlist = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    prefs = pd.read_csv('./results/prefs/unscaled_prefs.csv')
    ent_prefs = dms_tools2.prefs.prefsEntropy(prefs, charlist)
    ent_prefs.to_csv(entropy_file)

def MapDiffSeltoPDB(infile, 
                     scriptfile, 
                     colors = ['#fafafa', '#ff0000'], 
                     map_type = 'site_diffsel', 
                     restrict_to_chain = False, 
                     script_preamble = None,
                     script_postamble = False,
                     condition = None):
    '''Writes a colormapping script to be run in pymol; the colormapping is based on fracsurvive 
    to color a structure'''
    df = pd.read_csv(infile)
    df = df.dropna()
    column_names = list(df)
    #print(df)
        
    #print(df)
    
    # establish the color spectrum in hex and rgb.
    n_subdivisions = 500 # the color spectrum will be divided into this many discrete colors
    color1 = Color(colors[0])
    color2 = Color(colors[1])
    rgb_spectrum = [c.rgb for c in color1.range_to(color2, n_subdivisions)]
    rgb_spectrum_dict = dict([(i, rgb_spectrum[i]) for i in range(len(rgb_spectrum))])
    
    if map_type == 'site_diffsel':
        assert 'positive_diffsel' in column_names
        min_avg = df.min()['positive_diffsel']  
        max_avg = df.max()['positive_diffsel']  # the min and max will be mapped to color1 and color2, respectively
        range_avg = max_avg - min_avg
        df['colorindex'] =  (df.positive_diffsel - min_avg)/range_avg*(n_subdivisions-1)
        
    elif map_type == 'entropy':
        assert 'entropy' in column_names
        min_val = df.min()['entropy']  
        max_val = df.max()['entropy']  # the min and max will be mapped to color1 and color2, respectively
        range_val = max_val - min_val
        df['colorindex'] =  (df.entropy - min_val)/range_val*(n_subdivisions-1)
    
    df['colorindex'] = df['colorindex'].astype(int) # round to nearest index
    df['rgb'] = df['colorindex'].map(rgb_spectrum_dict)        
    site_color_mapping = pd.concat([df['site'], df['rgb']], axis=1)

    # write out the script to *scriptfile*:
    f = open(scriptfile, 'w')
    
    if script_preamble:
        preamblef = open(script_preamble, 'r')
        for line in preamblef.readlines():
            f.write(line)
        f.write('\n\n')
        preamblef.close()
    
    for i in range(len(df.index)):
        rgblist = [min(1, c) for c in site_color_mapping.iloc[i]['rgb']]
        r = site_color_mapping.iloc[i]['site']
        
        f.write("cmd.set_color(\'color{0}\', \'{1}\')\n".format(r, rgblist))
        f.write("cmd.color(\'color{0}\', \'resi {0}\')\n".format(r))
    
    if script_postamble:
        f.write("""cmd.create("dimer", "chain C+E")
cmd.show("surface", "dimer")
cmd.hide("everything", "5IRE")
cmd.set_view ("\\
    -0.352205604,    0.911581576,   -0.212054998,\\
     0.782968819,    0.162855446,   -0.600361347,\\
    -0.512743294,   -0.377484351,   -0.771101713,\\
     0.000000000,    0.000000000, -447.172607422,\\
  -115.573463440, -110.383399963, -125.525375366,\\
   352.554290771,  541.790893555,  -20.000000000" )

cmd.select("glycan1", "/dimer//C/NAG`600+601")
cmd.select("glycan2", "/dimer//E/NAG`600+601")
cmd.hide("everything", "glycan1")
cmd.hide("everything", "glycan2")
cmd.deselect()
cmd.set("ray_opaque_background", "off")
cmd.png("MR766_{0}_top.png", ray=1)

cmd.set_view ("\\
    -0.375719965,   -0.311803222,   -0.872701287,\\
     0.787538230,   -0.603787839,   -0.123329692,\\
    -0.488471270,   -0.733628869,    0.472416550,\\
     0.000000000,    0.000000000, -447.172607422,\\
  -115.573463440, -110.383399963, -125.525375366,\\
  -21650.402343750, 22544.750000000,  -20.000000000" )

cmd.png("MR766_{0}_side.png", ray=1)""".format(condition))
        f.write('\n\n')
    f.close()


MapDiffSeltoPDB("./results/diffsel/summary_ZKA64-meansitediffsel.csv", 
                 './pymol_scripts/ZKA64_mean_sitediffsel.py', 
                 map_type = 'site_diffsel',
                 script_preamble = False,
                 script_postamble = True,
                 condition = 'ZKA64')


MapDiffSeltoPDB("./results/diffsel/summary_ZKA185-meansitediffsel.csv", 
                 './pymol_scripts/ZKA185_mean_sitediffsel.py', 
                 map_type = 'site_diffsel',
                 script_preamble = False,
                 script_postamble = True,
                 condition = 'ZKA185') 

MapDiffSeltoPDB("./results/prefs/rescaled_prefs_entropy.csv", 
                 './pymol_scripts/MR766_entropy.py', 
                 map_type = 'entropy',
                 colors = ['#fafafa', '#007713'],
                 script_preamble = False,
                 script_postamble = True,
                 condition = 'entropy')


