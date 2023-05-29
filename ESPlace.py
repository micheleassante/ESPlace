import argparse
import pandas as pd
from solv_repo import *
from cubegen import * 
from footprinting_mod import footprinter
from explicit_solv_mod import *
# pd.options.mode.chained_assignment = None





parser = argparse.ArgumentParser(prog='ESPlace', description=""" 
ESPlace takes a gaussian log file as input, finds checkpoint file relative to that log file 
and generates an electron-density and an electrostatic potential cube file in ESP surfaces. From these two files, 
it identifies the most likely interaction points around the solutes. These can be optionally used to automatically
place explicit solvent molecules. \n 

Supported solvents: Water, DMSO, Dichlorobenzene, 1,4-Dioxane, Acetophenone, Anisole, Benzonitrile,
Butanone, Chloroform, Diethylether, EthylAcetate, THF, Methanol, Toluene""")
parser.add_argument('filename', help = 'input is a gaussian log file') 
parser.add_argument('-m','--surface', help="Print full ESP surface mesh", required=False, action='store_true')
parser.add_argument('-f','--footprint', help="Generate footprint of solute (solute + interaction points)", required=False, action='store_true')
parser.add_argument('-n','--number', help='Number of interaction points and thus explicit solvent molecule requested. format: (max, min)', default=(5,5))
parser.add_argument('-x','--explicit', help="Create input geometry with explicit solvent molecules", required=False, action='store_true')
parser.add_argument('-s','--solvent', help="Select type of solvent, default is water", default='water')
parser.add_argument('-i','--isodensity', help="Select range value to be used as filter for isodensity level", default=0.002)
parser.add_argument('-r','--radius', help="Select value for radius of footprint", default=1.74)
args = vars(parser.parse_args())


# read gaussian log file as input and extracts level of theory and name of check point file
input_log = args['filename']
name = input_log.split('.')[0]
chk_file, lev_theo = read_logfile(input_log)



# from checkpoint file calculate, fchk file, density and electrostatic potential cube
fchk_file = generate_fchk(chk_file)
den_file = generate_dencub(fchk_file)
esp_file = generate_espcub(fchk_file)


iso_filter = float(args['isodensity'])
requested_int = [int(args['number'].split(',')[0]),  int(args['number'].split(',')[1])]
radius = float(args['radius'])


solvent_dict = globals()[args['solvent']]
solvent_name = args['solvent']


# function to print gaussian input files
def gwriter(name, geometry, mem = '1GB', proc='1', levtheo='# m06/6-31+G** opt freq integral(grid=ultrafine)'):

    with open(name, 'w') as f:
        f.write('{}\n'.format('%mem=' + mem))
        f.write('{}\n'.format('%nproc=' + proc))
        f.write('{}\n'.format(r'%chk=' + name.strip('.com') + '.chk'))
        f.write('{}\n'.format(levtheo))
        
        f.write('{}\n'.format('Title'))
        f.write('\n')
        f.write('{}\n'.format('0 1'))
        for atom in geometry:
            s = str(atom).replace('[', '').replace(']', '').replace("'", '')
            f.write(s + '\n')
        f.write('\n')

# generate solute geometry, surface mesh, interaction points and footprint
solute, surface, i_points, footprint  = footprinter(den_file, esp_file, iso_filter, radius, requested_int )

# if requested prints surface to csv file
if args['surface'] == True:
    surface.to_csv(name + "_surface.csv", header = None, index = False)

# if requested prints footprint to gaussian input file for visualization
if args['footprint'] == True:
    gwriter(name + '_fp.com', footprint, '1GB', '1', lev_theo)

# if requested creates mixed implicit/explicit geometry and writes it to gaussian input file
if args['explicit'] == True:
    cluster = cluster_maker(solute, i_points, solvent_dict, solvent_name)
    gwriter(name + '_hybrid.com', cluster, '1GB', '1', lev_theo)

