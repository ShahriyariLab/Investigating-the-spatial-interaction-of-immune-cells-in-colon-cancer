'''
Investigating the spatial interaction of immune cells in colon cancer

Smoothen: Smoothens discountinous data that are projected onto a continuous space

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
def smoothen(Th0,Tc0,Tr0,Dn0,M0,C0,S1,mesh):

    ###############################################################
    #Function spaces are chosen so that they match the subspaces from the main script
    ###############################################################
    P1 = FiniteElement('P', triangle,1)
    PB = FiniteElement('B', triangle,3)
    NEE = FiniteElement('P', triangle,1)
    element = MixedElement([NEE,NEE,NEE,NEE,NEE,NEE])
    Mixed_Space = FunctionSpace(mesh, element)
    ###############################################################

    ###############################################################
    #Functions
    ###############################################################
    U = Function(Mixed_Space)
    Th,Tc,Tr,Dn,M,C= split(U)
    v1, v2, v3, v4, v5, v6= TestFunctions(Mixed_Space)
    ############################################################################

    ###############################################################
    #Fake diffusion coefficients
    ###############################################################
    D_1 = 5e-5
    ###############################################################

    ###############################################################
    #Diffusion equation weak form and solve
    ###############################################################
    F1 = D_1*dot(grad(Th),grad(v1))*dx+Th*v1*dx-Th0*v1*dx\
    + D_1*dot(grad(Tc),grad(v2))*dx+Tc*v2*dx-Tc0*v2*dx\
    + D_1*dot(grad(Tr),grad(v3))*dx+Tr*v3*dx-Tr0*v3*dx\
    + D_1*dot(grad(Dn),grad(v4))*dx+Dn*v4*dx-Dn0*v4*dx\
    + D_1*dot(grad(M),grad(v5))*dx+M*v5*dx-M0*v5*dx\
    + D_1*dot(grad(C),grad(v6))*dx+C*v6*dx-C0*v6*dx\

    solve(F1==0,U,[])
    ###############################################################

    ###############################################################
    #Splitting and projecting for return
    ###############################################################
    Th_,Tc_,Tr_,Dn_,M_,C_= U.split()
    Th = project(Th_,S1)
    Tc = project(Tc_,S1)
    Tr = project(Tr_,S1)
    Dn = project(Dn_,S1)
    M = project(M_,S1)
    C = project(C_,S1)
    ###############################################################
    return Th,Tc,Tr,Dn,M,C
###############################################################
