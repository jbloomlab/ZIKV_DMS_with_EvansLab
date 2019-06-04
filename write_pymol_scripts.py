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
                     script_postamble = True,
                     condition = True):
    '''Writes a colormapping script to be run in pymol to color a structure'''
    df = pd.read_csv(infile)
    df = df.dropna()
    column_names = list(df)

    
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
        postamblef = open(script_postamble, 'r')
        for line in postamblef:
            if "__" in line:
                line = line.replace("__", "_{0}_".format(condition))
            f.write(line)
        postamblef.close()




structures = ['5ire', '5u4w', '6co8']

for structure in structures: 
    postamble_file = 'pymol_scripts/'+structure+'_postamble.txt'

    MapDiffSeltoPDB("./results/diffsel/summary_ZKA64-meansitediffsel.csv", 
                    './pymol_scripts/{0}_ZKA64_mean_sitediffsel.py'.format(structure), 
                    map_type = 'site_diffsel',
                    script_preamble = False,
                    script_postamble = postamble_file,
                    condition = "ZKA64")


    MapDiffSeltoPDB("./results/diffsel/summary_ZKA185-meansitediffsel.csv", 
                    './pymol_scripts/{0}_ZKA185_mean_sitediffsel.py'.format(structure), 
                    map_type = 'site_diffsel',
                    script_preamble = False,
                    script_postamble = postamble_file,
                    condition = "ZKA185") 

    MapDiffSeltoPDB("./results/prefs/rescaled_prefs_entropy.csv", 
                    './pymol_scripts/{0}_MR766_entropy.py'.format(structure), 
                    map_type = 'entropy',
                    colors = ['#fafafa', '#007713'],
                    script_preamble = False,
                    script_postamble = postamble_file,
                    condition = "entropy")


