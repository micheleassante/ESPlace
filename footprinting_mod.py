
import pandas as pd
import copy 
from itertools import zip_longest



def footprinter(dencub, espcub, filter, radius, number_exp = [1,1]):
    
    aa=1.0/0.5291772083 # convert distances to anstroms  

    solute_coord=[]
    spacing_vec=[]
    nline=0
    values=[]
    esp=[]
    origin = []
    # Read den cube file and parse all data
    for line in open(dencub,"r"):
        nline+=1
        if nline==3:
            try:
                nat=int(line.split()[0])
                origin=[float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]
            except:
                print("ERROR: non recognized cube format")
        elif nline >3 and nline <= 6:
            spacing_vec.append(line.split())
        elif nline > 6 and nline <= 6+nat:
            solute_coord.append([line.split()[0], float(line.split()[2])/aa,float(line.split()[3])/aa,float(line.split()[4])/aa])
        elif nline > 5:
            if nline > 6+nat:
                for i in line.split():
                    values.append(float(i))

    nline=0
    for line in open(espcub,"r"):
        nline+=1
        if nline==3:
            try:
                cat=int(line.split()[0])
                origin=[float(line.split()[1])/aa,float(line.split()[2])/aa,float(line.split()[3])/aa]
            except:
                print("ERROR: non recognized cube format")

        elif nline > 6 and nline <= 6+nat:
            pass
        elif nline > 5:
            if nline > 6+cat:
                for i in line.split():
                    esp.append(float(i))
        else:
            pass

    idx=-1

    #store x,y,z,den,esp values in data 
    data=[]
    for i in range(0,int(spacing_vec[0][0])):
        for j in range(0,int(spacing_vec[1][0])):
            for k in range(0,int(spacing_vec[2][0])):
                idx+=1
                x,y,z= i*float(spacing_vec[0][1]),j*float(spacing_vec[1][2]),k*float(spacing_vec[2][3])
                data.append([x/aa,y/aa,z/aa, values[idx], esp[idx]])


    # values into a 4 columns array in pandas with labeled columns
    df = pd.DataFrame(data, columns=['x','y','z','Density','ESP'])
    
    # center the dataframe
    df['xc'] = df['x'] + origin[0]
    df['yc'] = df['y'] + origin[1]
    df['zc'] = df['z'] + origin[2]

    n = filter
    #create data frame with values of Density between 90 and 110% of given isodensity value
    vdw = df[(df['Density'] >= (n/100)*90) & (df['Density'] <= (n/100)*110)]
    coord  = vdw[['xc','yc','zc', 'ESP']]
    surface_esp = copy.deepcopy(coord)
    # find ESP interaction points 
    r = radius # value from Hunter paper
    ssips = pd.DataFrame()

    max_point_range = number_exp[0]
    min_point_range = number_exp[1]
    
    # iterate through the selected number of MAX and MIN points
    for max_pt, min_pt in zip_longest(range(max_point_range), range(min_point_range)):
        
        if max_pt is not None:
            try:
                #find maximum 
                top1 = coord.nlargest(1, 'ESP')
                # carve the surface around the maxima and store the coordinates of maxima + value                    
                c1,c2,c3 = float(top1['xc']), float(top1['yc']), float(top1['zc'])                    
                ssips=ssips.append(top1)
                coord = coord.loc[((coord["xc"] - c1)**2) + ((coord["yc"] - c2)**2) + ((coord["zc"] - c3)**2) > (r**2)]
            
            except:
                pass
    
        if min_pt is not None:
            try:
                #find minimum
                bot1 = coord.nsmallest(1, 'ESP')
                # carve the surface around the maxima and store the coordinates of maxima + value      
                b1,b2,b3 = float(bot1['xc']), float(bot1['yc']), float(bot1['zc'])
                ssips=ssips.append(bot1)
                coord = coord.loc[((coord["xc"] - b1)**2) + ((coord["yc"] - b2)**2) + ((coord["zc"] - b3)**2) > (r**2)]
            
            except:
                pass

    # save interaction points in dictionary with ESP as keys and coordinates as values
    int_point2 = ssips.values.tolist()
    
    i_pts = copy.deepcopy(int_point2)

    # create footprint geometry
    footprint=[]
    
    for i in solute_coord:
        footprint.append(i)
    for j in int_point2:
        if j[-1] > 0:
            j.insert(0, 'He')
            footprint.append(j)
        else:
            j.insert(0, 'X')
            footprint.append(j)
    
    return solute_coord, surface_esp, i_pts, footprint  
