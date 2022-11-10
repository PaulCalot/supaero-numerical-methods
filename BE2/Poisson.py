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
h = 0.04 # to change the default resolution of the mesh (default defined in Flag.geo file)
os.system('gmsh -2 Flag.geo -setnumber h '+str(h))

m = gf.Mesh('import','gmsh','Flag.msh')

print('')
print('===')

print(m.regions(),'region ids from .geo file (not the Physical ones!)')
print(m.memsize(),'bytes used to store the mesh')

m.display() # To get some infos on the mesh

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
mfu = gf.MeshFem(m,1) # The FE space is on mesh m, for a scalar field
mfu.set_fem(gf.Fem('FEM_PK(2,2)')) # PK-Lagrange, continuous, on simplex of dim. 2, order 2
# mfu.set_classical_fem(2) # Do = PK-Lagrange, continuous, on simplex of dim. 2, order 2

# print(help(mfu)) # This object represents a finite element method defined on a whole mesh.

# An exact integration associated to the mesh
mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(4)')) # integration pol. of order 4 is exact

# print(help(mim)) # This object represents an integration method defined on a whole mesh (an 
#  |  potentially on its boundaries).



#%% Model definition

# Data (may vary in x and y)
# .eval : interpolate an expression on the (lagrangian) MeshFem.
f = mfu.eval('-1.') # RHS 
gn = mfu.eval('0.') # Neumann
gd_l = mfu.eval('0.') # Dirichlet
gd_r = mfu.eval('-1.') # Dirichlet

# Model, variables
md = gf.Model('real') # Can handle complex fields
# print(help(md)) #     Model variables store the variables and the state data and the
#  |  description of a model. This includes the global tangent matrix, the right
#  |  hand side and the constraints. There are two kinds of models, the `real`
#  |  and the `complex` models.
md.add_fem_variable('u', mfu) # Name the variable to compute

# Data initialization
md.add_initialized_fem_data('f', mfu, f)
md.add_initialized_fem_data('gn', mfu, gn)
md.add_initialized_fem_data('gd_l', mfu, gd_l)
md.add_initialized_fem_data('gd_r', mfu, gd_r)

#%% Bricks

#md.add_Laplacian_brick(mim, 'u') # Laplacian brick
md.add_linear_term(mim, 'Grad_u:Grad_Test_u') # Laplacian variational formulation, produit contracté
# c'est la fonction a je pense

md.add_source_term_brick(mim, 'u', 'f') # RHS terms in the variational formulation ?
md.add_source_term_brick(mim, 'u', 'gn', BOTTOM_BOUND) # on ajoute les termes sources uniquement sur les bords 
                                                       # (i.e. à quoi doit être égal u sur ces domaines, ici 1D)
md.add_source_term_brick(mim, 'u', 'gn', TOP_BOUND)
md.add_source_term_brick(mim, 'u', 'gn', HOLE_BOUND)

# Essential conditions with simplification (other ways for Dirichlet imposition exist)
md.add_Dirichlet_condition_with_simplification('u', LEFT_BOUND, 'gd_l')
md.add_Dirichlet_condition_with_simplification('u', RIGHT_BOUND, 'gd_r')

md.brick_list() # Get some infos on the model now defined

#%% Solve with the simplest way

md.solve()
u = md.variable('u') # extraction of the computed solution
mfu.export_to_pos('Laplacian.pos',mfu,u,'u') # save as gmsh .pos file for visualization

#%% Solve with algebra

A = md.tangent_matrix()
A.to_csc()
b = md.rhs()

u_mumps = np.transpose(gf.linsolve('mumps',A,b))
mfu.export_to_pos('Laplacian_mumps.pos',mfu,u_mumps,'u (mumps)')

error = u-u_mumps
print('L-sup error w.r.t GetFEM solution: ',np.max(np.abs(error)))
