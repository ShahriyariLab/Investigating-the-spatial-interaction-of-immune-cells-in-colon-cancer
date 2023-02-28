'''
Investigating the spatial interaction of immune cells in colon cancer

IC_Loc_DG: Creating spatially distribute initial conditions based on biological
            data

Author: Navid Mohammad Mirzaei https://sites.google.com/view/nmirzaei
                               https://github.com/nmirzaei

(c) Shahriyari Lab https://sites.google.com/site/leilishahriyari/
'''
###############################################################
#Importing required functions
###############################################################
from fenics import *
import numpy as np
import pandas as pd
import csv
import random as rnd
###############################################################
def IC_Loc_DG(mesh,S1):

    ###############################################################
    #Reads the spatial data
    ###############################################################
    IC_Th = pd.read_csv('input/input/Th.csv').to_numpy()
    IC_Tc = pd.read_csv('input/input/Tc.csv').to_numpy()
    IC_Tr = pd.read_csv('input/input/Tr.csv').to_numpy()
    IC_Dn = pd.read_csv('input/input/Dn.csv').to_numpy()
    IC_M = pd.read_csv('input/input/M.csv').to_numpy()
    IC_C = pd.read_csv('input/input/C.csv').to_numpy()
    ###############################################################
    #Creating bounding box tree to be able to check what cell belongs to which mesh triangle
    ###############################################################
    X = mesh.coordinates()
    tree = BoundingBoxTree()
    tree.build(mesh)
    Num_cell = mesh.num_cells()
    ############################################################################

    ###############################################################
    #Functions on the fine mesh
    ###############################################################
    Th = Function(S1)
    Tc = Function(S1)
    Tr = Function(S1)
    Dn = Function(S1)
    M  = Function(S1)
    C  = Function(S1)
    ############################################################################

    ###############################################################
    #Counting cells in each triangle based on the spatial data
    ###############################################################
    id = 10000000000
    for i in range(len(IC_Th)):
        p = Point(IC_Th[i,0]-0.3,IC_Th[i,1]-0.225,0)
        id = tree.compute_first_entity_collision(p)
        if id>Num_cell:
            continue
        else:
            Th.vector()[id] = Th.vector().get_local()[id]+1

    id = 10000000000
    for i in range(len(IC_Tc)):
        p = Point(IC_Tc[i,0]-0.3,IC_Th[i,1]-0.225,0)
        id = tree.compute_first_entity_collision(p)
        if id>Num_cell:
            continue
        else:
            Tc.vector()[id] = Tc.vector().get_local()[id]+1

    id = 10000000000
    for i in range(len(IC_Tr)):
        p = Point(IC_Tr[i,0]-0.3,IC_Tr[i,1]-0.225,0)
        id = tree.compute_first_entity_collision(p)
        if id>Num_cell:
            continue
        else:
            Tr.vector()[id] = Tr.vector().get_local()[id]+1

    id = 10000000000
    for i in range(len(IC_Dn)):
        p = Point(IC_Dn[i,0]-0.3,IC_Dn[i,1]-0.225,0)
        id = tree.compute_first_entity_collision(p)
        if id>Num_cell:
            continue
        else:
            Dn.vector()[id] = Dn.vector().get_local()[id]+1

    id = 10000000000
    for i in range(len(IC_M)):
        p = Point(IC_M[i,0]-0.3,IC_M[i,1]-0.225,0)
        id = tree.compute_first_entity_collision(p)
        if id>Num_cell:
            continue
        else:
            M.vector()[id] = M.vector().get_local()[id]+1

    id = 10000000000
    for i in range(len(IC_C)):
        p = Point(IC_C[i,0]-0.3,IC_C[i,1]-0.225,0)
        id = tree.compute_first_entity_collision(p)
        if id>Num_cell:
            continue
        else:
            C.vector()[id] = C.vector().get_local()[id]+1

        ############################################################################
    ###############################################################

    ###############################################################
    #Max values
    ###############################################################
    maxTh = max(Th.vector()[:])
    maxTc = max(Tc.vector()[:])
    maxTr = max(Tr.vector()[:])
    maxDn = max(Dn.vector()[:])
    maxM = max(M.vector()[:])
    maxC = max(C.vector()[:])


    Th.vector()[:]/=maxTh
    Tc.vector()[:]/=maxTc
    Tr.vector()[:]/=maxTr
    Dn.vector()[:]/=maxDn
    M.vector()[:]/=maxM
    C.vector()[:]/=maxC
    ############################################################################

    ###############################################################
    #Save the IC plots
    ###############################################################
    File('DG_IC/Th.pvd')<<Th
    File('DG_IC/Tc.pvd')<<Tc
    File('DG_IC/Tr.pvd')<<Tr
    File('DG_IC/Dn.pvd')<<Dn
    File('DG_IC/M.pvd')<<M
    File('DG_IC/C.pvd')<<C
    ############################################################################


    return Th,Tc,Tr,Dn,M,C,maxTh,maxTc,maxTr,maxDn,maxM,maxC
###############################################################
