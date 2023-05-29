import numpy as np
from numpy.random.mtrand import f
from scipy.spatial.transform import Rotation as R


# function for the identification of closest center of atoms within a 7Å radius
def find_closest_centroid(mol, i_point, rad=7):
    neigh = []
    for e in mol:
        rf = np.linalg.norm(e-i_point)
        if rf < rad:
            neigh.append([e[0],e[1],e[2]])
    
    neighbours = np.array([np.array(x) for x in neigh]) 
    close_centroid = np.array((sum(neighbours[:,0])/len(neighbours[:,0]),sum(neighbours[:,1])/len(neighbours[:,1]),sum(neighbours[:,2])/len(neighbours[:,2])))
    
    return close_centroid


def solvent_placer(solv_repo, ssip, solute_cent):
    keys = [l[0] for l in solv_repo]
    solvent_array = np.array([l[1] for l in solv_repo])

    vo = solvent_array[-2]
    xo = ssip
    solvent_o = np.array([at-vo for at in solvent_array]) # translate standard solvent coordinates at the origin


    
  
    
    ox, oy, oz = solute_cent 
    xx, xy, xz = xo - solute_cent
    Vx, Vy, Vz = solvent_o[-2]
    bx, by, bz = solvent_o[-1]
    r = np.array([xx, xy, xz])
    s = np.array([bx-Vx, by-Vy, bz-Vz])
    rs = np.cross(r,s)
    
   



    rs2 = rs / np.sqrt(np.dot(rs, rs))
    
    alfa = np.arccos((np.dot(r,s)/(np.linalg.norm(r) * np.linalg.norm(s))))

    # translation + first rotation
    
    r1 = R.from_rotvec(alfa * rs2)
    rot1 = r1.apply(solvent_o)

    if np.allclose(r, rot1[-1], atol=1e-01) == False:
        r_neg = R.from_rotvec(-alfa * rs2)
        rot1 = r_neg.apply(solvent_o)
    else:
        pass
    
    betas = np.arange(0,3.14, 0.5)
    
    rotamer = []
    distance = 0
    for b in betas:
        ac = rot1[-1] - rot1[-2]
        l = ssip - rot1[-2]
        r2 = R.from_rotvec(b * ac)
        rot2 = r2.apply(rot1)

        
        finalfinal = np.array([at + l for at in rot2])


        
        fdist = [np.linalg.norm(k - solute_cent) for k in finalfinal[0:-2]]
        avg_dist = sum(fdist)/len(fdist)
        
        if avg_dist > distance:
            distance = avg_dist
            rotamer = finalfinal
        else:
            pass
    
    ex1 = [[str(t) for t in l] for l in rotamer]
    for i in range(len(ex1)):
        rot_atom = ex1[i]
        rot_atom.insert(0, keys[i])
         
    
    
 
    solv = [e for e in ex1 if e[0] not in ('Cnt', 'Vx')]
    return solv


def cluster_maker(solute, int_points, sv_dict, sv_name):
    hybrid_solute = [[i for i in d] for d in solute]    
    for key in int_points:
        esp = key[3] # esp value of interaction points 
        ssip = np.array(key[0:3]) # coordinates of interaction points
        
        #create numpy array of solute coordinates and localize the centeroid for optimal solvent orientation
        solute_coord = [[float(c) for c in i[1:]] for i in hybrid_solute]
        at_array = np.array([np.array(x) for x in solute_coord]) 
        closest = find_closest_centroid(at_array, ssip) # finds centroid in a 5A radius from the inter. point
    
        # place solvent molecules according to ESP value of the interaction point and centroid position. 
        # Newly placed solvent molecule is added to solute atoms array and its atoms are considered in
        # the identification of the new centroid for the successive solvent placement iteration
        if esp > 0.0001:
            solvent_mode = sv_dict[sv_name + '_min']
            exp_solv = solvent_placer(solvent_mode, ssip, closest)
            for q in exp_solv:
                hybrid_solute.append(q)
        

        elif esp < -0.0001:
            solvent_mode = sv_dict[sv_name + '_max']
            exp_solv = solvent_placer(solvent_mode, ssip, closest)
            for q in exp_solv:
                hybrid_solute.append(q)
            
        else:
            pass

    hybrid_solute
    
    return hybrid_solute