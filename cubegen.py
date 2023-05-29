
import subprocess
import re
import os



# read gaussian log file and return name of chk and level of theory
def read_logfile(logfile):
    with open(logfile, 'r') as f:
        pattern_chk = r'%chk'
        pattern_lev = '#'
        lines = f.readlines()
        
        chk_file = ''
        lev_theo = ''

        for line in lines:
            try:
                if re.search(pattern_chk, line):  # Check if the line contains the pattern
                    chk_file = line.split('=')[1]
                if re.search(pattern_lev, line):
                    lev_theo = line.lstrip()
                    break
            except:
                print('Cannot find checkpoint file. Please make sure a check point file is available')
        
    
        
        return chk_file, lev_theo

### run formchl to get formatted checkpoint 

def generate_fchk(input_chk):
    pwd = os.getcwd()
    fchk_out = input_chk.replace('chk', 'fchk')
    fchk_path = os.path.join(pwd, fchk_out)

    if os.path.exists(fchk_path):
        return fchk_out 
    
    else:
        cmd_fchk= 'formchk %s' % (input_chk)
        subprocess.run(cmd_fchk.split())
        return fchk_out 


### run cubegen for density

def generate_dencub(input_name):
    
    pwd = os.getcwd()
    clean_name = input_name.split('.')[0]
    den_out = '%s_den.cub' % (clean_name) 
    den_path = os.path.join(pwd, den_out)    
    
    if os.path.exists(den_path):
        return den_out
    else:
        print('... Calculating electron density cube ...')
        cmd_den = 'cubegen 0 fdensity=scf %s %s_den.cub' % (clean_name, clean_name)
        subprocess.run(cmd_den.split())
        den_out = '%s_den.cub' % (clean_name)
        return den_out

## run cubegen for potential
def generate_espcub(input_name):
    
    pwd = os.getcwd()
    clean_name = input_name.split('.')[0]
    esp_out = '%s_esp.cub' % (clean_name) 
    esp_path = os.path.join(pwd, esp_out)   
    
    if os.path.exists(esp_path):
        return esp_out
    
    else:
        print('... Calculating electrostatic potential cube ...')
        cmd_esp = 'cubegen 0 potential=scf %s %s_esp.cub' % (clean_name, clean_name)
        subprocess.run(cmd_esp.split())
        esp_out = '%s_esp.cub' % (clean_name)
        return esp_out

