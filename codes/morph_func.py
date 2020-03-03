from allensdk.api.queries.biophysical_api import BiophysicalApi
from allensdk.api.queries.cell_types_api import CellTypesApi

import os.path

from pandas import Series, DataFrame
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import os
import h5py
from scipy import signal
import pylab as pl
from sklearn.decomposition import PCA
import numpy.linalg as la
import math
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from allensdk.core.cell_types_cache import CellTypesCache
import matplotlib as mpl
from mpl_toolkits import mplot3d

# get the rotation angle for specific cellid
def get_rotation_theta(cell_id):
    # get the morphology data (from apical dendrite, dendrite, and soma) used for PCA
    [morph_data,morph_soma]=get_cell_morphXYZ(cell_id)
    [v,theta] = cal_rotation_angle(morph_data)
    return theta


# Obtain morphlogoy location for PCA     
def get_cell_morphXYZ(cell_id):
    
    from allensdk.core.cell_types_cache import CellTypesCache
    from allensdk.core.swc import Marker
    import pprint

    ctc = CellTypesCache(manifest_file='cell_types/manifest.json')
    morph = ctc.get_reconstruction(cell_id) 
    
    markers = ctc.get_reconstruction_markers(cell_id) 
 
    x=[]
    y=[]
    z=[]  
    for n in morph.compartment_list:
        #print(n['type']) #type=1, soma; type=2, axon; type=3, dendrite; type=4,apical dendrite 
        if n['type']==4 or n['type']==3 or n['type']==1:
            x.append(n['x']-morph.soma["x"])
            y.append(n['y']-morph.soma["y"])
            z.append(n['z']-morph.soma["z"])
    
    morph_data = np.array(np.column_stack((x,y,z))) 
    morph_soma = [morph.soma["x"],morph.soma["y"],morph.soma["z"]]
    
    return (morph_data,morph_soma)


# RXYZ, the same as the bmtk
def cal_rotation_angle(morph_data):
    
    pca = PCA(n_components=2)
    pca.fit(morph_data)
    proj = morph_data.dot(pca.components_[0])  # the projection of morphology on the direction of first pca
    
    #v1 = -1*pca.components_[0]  # the first principal component, when apical dendrite goes down 
    #v1 = 1*pca.components_[0]  # the first principal component 
    
    v1 = np.sign(proj.mean())*pca.components_[0]
    
    # The goal is to rotate v1 to parallel to y axis
    x1=v1[0]
    y1=v1[1]
    z1=v1[2]
    
    # First rotate in the anticlockwise direction around z axis untill x=0
    v2 = [0,math.sqrt(y1*y1+x1*x1),z1]
    dv = [v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]]
    anglez= 2*math.asin(math.sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2])*0.5/v2[1])
    if x1<0:  # when x1 in the negative side, change the sign
        anglez=-anglez

    # Second rotate in the anticlockwise direction round x axis untill z = 0
    v3 = [0,math.sqrt(x1*x1+y1*y1+z1*z1),0]
    dv2= [v3[0]-v2[0],v3[1]-v2[1],v3[2]-v2[2]]
    anglex= -2*math.asin(math.sqrt(dv2[0]*dv2[0]+dv2[1]*dv2[1]+dv2[2]*dv2[2])*0.5/v3[1])
    if z1<0:  # when z1 in the negative side, change the sign
        anglex=-anglex

    theta=[anglex,0,anglez]
    R = eulerAnglesToRotationMatrix(theta)
    v3_hat=R.dot(v1)

    return (v1,theta) 

# Calculates Rotation Matrix given euler angles.
def eulerAnglesToRotationMatrix(theta):

    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                    ])
 
    R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                    ])
                 
    R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])
    R = np.dot(R_x, np.dot( R_y, R_z ))

    return R

# function to plot morphology
def cell_morphology_rot(cell_id, x_soma, y_soma, z_soma, theta):

    theta_z = theta[2]
    theta_y = theta[1]
    theta_x = theta[0]
    
    # download and open an SWC file
    ctc = CellTypesCache(manifest_file='cell_types/manifest.json')
    morph = ctc.get_reconstruction(cell_id) 
    
    #First applying a rotation angle around z axis
    tr_rot_z = [math.cos(theta_z),-math.sin(theta_z),0,
                math.sin(theta_z),math.cos(theta_z),0,
                0,0,1,
                0,0,0
              ]   
    
    #Second applying a rotation angle around y axis
    tr_rot_y = [math.cos(theta_y),0,math.sin(theta_y),
                0,1,0,
                -math.sin(theta_y),0,math.cos(theta_y),
                0,0,0
              ]
    
    #Third applying a rotation angle around x axis
    tr_rot_x = [1,0,0,
                0,math.cos(theta_x),-math.sin(theta_x),            
                0,math.sin(theta_x),math.cos(theta_x),
                0,0,0
              ]

    
    morph.apply_affine(tr_rot_z)
    morph.apply_affine(tr_rot_y)
    morph.apply_affine(tr_rot_x)
       
    # translate the soma location
    tr_soma = [1, 0, 0,
               0, 1, 0,
               0, 0, 1,
               -morph.soma["x"]+x_soma, -morph.soma["y"]+y_soma, -morph.soma["z"]+z_soma
             ]
    morph.apply_affine(tr_soma)
        
    # plot xy and xz views
    #   fig, axes = plt.subplots(1, 2, sharey=True, sharex=True)
    #    vm.plot_morph_swc(morph,axes,color='r')
    
   # fig, axes = plt.subplots(1, 2, sharey=True, sharex=True)

    # Make a line drawing of x-y and y-z views    
    return morph


def plot_cell_morph_xyzy(axes,morph):
    soma_col=[134.0/255.0,134.0/255.0,148.0/255.0]
    axon_col=[93.0/255.0,127.0/255.0,177.0/255.0]
    dend_col=[153.0/255.0,40.0/255.0,39.0/255.0]
    apical_dend_col=[227.0/255.0,126.0/255.0,39.0/255.0]
    ap = 1
    
    for n in morph.compartment_list:
        for c in morph.children_of(n):
            if n['type']==2:
                axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color=axon_col,alpha=ap)
                axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color=axon_col,alpha=ap)
            if n['type']==3:
                axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color=dend_col,alpha=ap)
                axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color=dend_col,alpha=ap)
            if n['type']==4:
                axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color=apical_dend_col,alpha=ap)
                axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color=apical_dend_col,alpha=ap)
            if n['type']==1: #soma
                axes[0].scatter(n['x'],n['y'],s=math.pi*(n['radius']**2),color=soma_col)
                axes[1].scatter(n['z'],n['y'],s=math.pi*(n['radius']**2),color=soma_col)   

    axes[0].set_ylabel('y')
    axes[0].set_xlabel('x')
    axes[1].set_xlabel('z')
    simpleaxis(axes[0])
    simpleaxis(axes[1])


def simpleaxis(ax):
    #Hide the right and top spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

   
def plot_mor(axes,cell_id,add_angle):

    # get the morphology data (from apical dendrite, dendrite, and soma) used for PCA
    [morph_data,morph_soma]=get_cell_morphXYZ(cell_id)
    print(morph_data.shape)

    # get the first principal vector and the rotation angle
    [v,theta] = cal_rotation_angle(morph_data)

    # Based on the rotation angle to calculate the rotation matrix R
    R = eulerAnglesToRotationMatrix(theta)  # rotation matrix

    # the first principal component before and after rotated   
    v = v*400
    v_rot = R.dot(v)   
    # The morphology locations used for PCA after rotations
    X_rot = np.array(morph_data)    # The rotated position of new x,y,z
    for i in range(0,len(X_rot)):
        X_rot[i,:]=R.dot(morph_data[i,:])
    
    print(cell_id)
    print(theta)
    
    # The location of soma, defined by the user
    x_soma=0
    y_soma=0
    z_soma=0
    # The original morphology before rotations
    theta0=[0,0,0]
    morph0 = cell_morphology_rot(cell_id,x_soma,y_soma,z_soma,theta0)
    
    # The morphology after rotations
    theta[2] = theta[2]+add_angle
    morph_rot = cell_morphology_rot(cell_id,x_soma,y_soma,z_soma,theta)

    plot_cell_morph_xy(axes,morph_rot)

