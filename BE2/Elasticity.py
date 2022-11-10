#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 15:24:47 2022

@author: ghaine
"""

#%% Initialization

import os
import numpy as np

# import basic modules
import getfem as gf

#%% Import gmsh file

# External call to gmsh client (a python API exists)
# with (re-)definition of the mesh size (see .geo file)
h = 0.04
os.system('gmsh -2 Flag.geo -setnumber h '+str(h))

m = gf.Mesh('import','gmsh','Flag.msh')

print('')
print('===')

print(m.regions(),'region ids from .geo file (not the Physical ones!)')
print(m.memsize(),'bytes used to store the mesh')

m.display()

#%% Prepare boundaries

# Merge the 4 regions for the hole
m.region_merge(11,12)
m.region_merge(11,13)
m.region_merge(11,14)
# Delete the 3 regions now in 11
m.delete_region([12,13,14])
print(m.regions(),'remaining region ids DIFFERENT from .geo file (definition of the "hole" boundary)')
# Name the indices for the sake of readability
BOTTOM_BOUND=15; RIGHT_BOUND=16; TOP_BOUND=17; LEFT_BOUND=18; HOLE_BOUND=11;
print('Bottom id: ', BOTTOM_BOUND, ', Right id: ', RIGHT_BOUND, ', Top id: ', TOP_BOUND,\
      ', Left id: ', LEFT_BOUND, ', Hole id: ', HOLE_BOUND)
# To be complete, we also name the domain
OMEGA=31;
print('Omega id: ', OMEGA)

#%% Finite element and integration method

# For the field u
mfu = gf.MeshFem(m,2) # The FE space is on mesh m, for a vector field of size 2
mfu.set_fem(gf.Fem('FEM_PK(2,2)')) # PK-Lagrange, continuous, on simplex of dim. 2, order 2

# An exact integration associated to the mesh
mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(4)')) # pol. of order 4 are exactly integrated

#%% Model definition

# Data
f = mfu.eval('[0., 0.]') # RHS
g1 = mfu.eval('[0., 0.]') # Neumann
g2 = mfu.eval('[1., 0.]') # Neumann
h = mfu.eval('[0., 0.]') # Dirichlet

E = 1e3
Nu = 0.3
Lambda = E*Nu/((1+Nu)*(1-2*Nu))
Mu = E/(2*(1+Nu))

# Model, variables, data
md = gf.Model('real') # Can handle complex fields
md.add_fem_variable('u', mfu)
md.add_initialized_data('cmu', Mu)
md.add_initialized_data('clambda', Lambda)
md.add_initialized_fem_data('f', mfu, f)
md.add_initialized_fem_data('g1', mfu, g1)
md.add_initialized_fem_data('g2', mfu, g2)
md.add_initialized_fem_data('h', mfu, h)

#%% Bricks


### TO DO ###


#%% Solve with the simplest way

md.solve()
u = md.variable('u')

mfdu=gf.MeshFem(m,1)
mfdu.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(2,1)'))
VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u','clambda','cmu', mfdu)

mfu.export_to_pos('Elasticity.pos',mfu,u,'u',mfdu,VM,'Von Mises')

#%% Solve with algebra 

A = md.tangent_matrix()
A.to_csc()
b = md.rhs()

u_mumps = np.transpose(gf.linsolve('mumps',A,b))
mfu.export_to_pos('Elasticity_mumps.pos',mfu,u_mumps,'u (mumps)')

error = u-u_mumps
print('L-sup error w.r.t GetFEM solution: ',np.max(np.abs(error)))
