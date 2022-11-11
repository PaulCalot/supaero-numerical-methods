#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 11:26:44 2022

@author: ghaine
"""

#  Initialization

import os
import numpy as np
import time
# import basic modules
import getfem as gf

#  Import gmsh file

# External call to gmsh client (a python API exists)
# with (re-)definition of the mesh size (see .geo file)
def run(epsilon=0.0001, Gamma = 20., h = 0.1): 
    os.system('gmsh -2 Flag.geo -setnumber h '+str(h))

    m = gf.Mesh('import','gmsh','Flag.msh')

    #print('')
    #print('===')

    #print(m.regions(),'region ids from .geo file (not the Physical ones!)')
    #print(m.memsize(),'bytes used to store the mesh')

    m.display()

    #  Prepare boundaries

    # Merge the 4 regions for the hole
    m.region_merge(11,12)
    m.region_merge(11,13)
    m.region_merge(11,14)
    # Delete the 3 regions now in 11
    m.delete_region([12,13,14])
    #print(m.regions(),'remaining region ids DIFFERENT from .geo file (definition of the "hole" boundary)')
    # Name the indices for the sake of readability
    BOTTOM_BOUND=15; RIGHT_BOUND=16; TOP_BOUND=17; LEFT_BOUND=18; HOLE_BOUND=11;
    #print('Bottom id: ', BOTTOM_BOUND, ', Right id: ', RIGHT_BOUND, ', Top id: ', TOP_BOUND,\
    #    ', Left id: ', LEFT_BOUND, ', Hole id: ', HOLE_BOUND)
    # To be complete, we also name the domain
    OMEGA=31;
    #print('Omega id: ', OMEGA)

    #  Finite element and integration method

    # For the field u
    mfu = gf.MeshFem(m,2) # The FE space is on mesh m, for a vector field of size 2
    mfu.set_classical_fem(2)

    mfp = gf.MeshFem(m,1) # The FE space is on mesh m, for a scalar field
    mfp.set_classical_fem(1)

    # An exact integration associated to the mesh
    mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(4)')) # pol. of order 4 are exactly integrated

    #  Model definition

    f = mfu.eval('[0., 0.]') # RHS
    h_inlet = mfu.eval('[4.*y*(1.-y), 0.]') # Dirichlet
    h_wall = mfu.eval('[0., 0.]') # Dirichlet

    nu = 0.001

    md = gf.Model('real')

    md.add_fem_variable('u', mfu)
    md.add_fem_variable('p', mfp)

    md.add_initialized_fem_data('f', mfu, f)
    md.add_initialized_fem_data('h_inlet', mfu, h_inlet)
    md.add_initialized_fem_data('h_wall', mfu, h_wall)
    md.add_initialized_data('nu', [nu])


    #  Bricks
    md.add_linear_term(mim,'nu * Grad_u:Grad_Test_u - p.Div_Test_u')
    md.add_linear_term(mim, 'Div_u.Test_p')
    md.add_linear_term(mim, 'nu * (Grad_u * Normal - u.Normal * Normal).Test_u', RIGHT_BOUND)
    # The following term could (should ?) work, however the pressures becomes very high
    # like it was not setting correctly the constant ?
    # md.add_linear_term(mim, 'p*Test_u.Normal - (nu * u.Normal * Test_u.Normal)', RIGHT_BOUND)

    # dirichlet
    md.add_source_term_brick(mim, 'u', 'f')
    md.add_Dirichlet_condition_with_simplification('u', LEFT_BOUND, 'h_inlet')
    md.add_Dirichlet_condition_with_simplification('u', HOLE_BOUND, 'h_wall') 
    md.add_Dirichlet_condition_with_simplification('u', BOTTOM_BOUND, 'h_wall') 
    md.add_Dirichlet_condition_with_simplification('u', TOP_BOUND, 'h_wall') 

    # md.add_nonlinear_term(mim, 'u.Grad_u.Test_u') # For incompressible Navier-Stokes

    #  Solve with the simplest way

    md.solve()

    u = md.variable('u');
    p = md.variable('p');

    mfu.export_to_pos('Stokes.pos',mfu,u,'u',mfp,p,'p')

    #  Solve with algebra 

    I_u = md.interval_of_variable('u') # Interval of dofs for u: start, nbdof
    Id_u = range(I_u[0],I_u[0]+I_u[1]) # Indices of dofs for u
    I_p = md.interval_of_variable('p') # Interval of dofs for p: start, nbdof
    Id_p = range(I_p[0],I_p[0]+I_p[1]) # Indices of dofs for p

    J = md.tangent_matrix()
    J.to_csc()
    b = md.rhs()

    A = gf.Spmat('copy',J,Id_u,Id_u)
    Bt = gf.Spmat('copy',J,Id_u,Id_p) # Warning, not exactly B^T as in the course because of the boundary condition on RIGHT_BOUND
    B = gf.Spmat('copy',J,Id_p,Id_u)
    Z = gf.Spmat('copy',J,Id_p,Id_p)

    F = b[Id_u]
    G = b[Id_p]

    # 
    # We use 'superlu' to solve our linear systems as it also
    # provides an estimate of the condition number
    # 

    #  Global Method

    t = time.time()

    X = gf.linsolve('superlu',J,b)
    U = np.transpose(X[0][Id_u]) # Warning [0] is needed (vector = first column of a matrix Nx1)
    P = np.transpose(X[0][Id_p]) # Warning [0] is needed (vector = first column of a matrix Nx1)

    elapsed = time.time() - t

    mfu.export_to_pos('Stokes_superlu.pos',mfu,U,'u (superlu)',mfp,P,'P (superlu)')

    error_u = u-U
    error_p = p-P
    #print('')
    #print('=== Global Method ===')
    #print('L-2 errors w.r.t GetFEM solution (Global): ',np.linalg.norm(error_u,2),' / ',np.linalg.norm(error_p,2))
    #print('L-2 norm of the divergence (Global): ',np.linalg.norm(B*U,2))
    #print('Estimation of the condition number (Global matrix): ', X[1])
    #print('Estimation of the memory usage (Global matrix): ', 8*8*3*J.nnz()*0.000008, 'Mb')
    #print('CPU Time: ', elapsed, ' s')

    tracking_global = {
        "method": "global",
        "duration": elapsed,
        "condition": X[1],
        "memory_usage": 8*8*3*J.nnz()*0.000008,
        "divergence_vector_field": np.linalg.norm(B*U,2),
        "error_u": np.linalg.norm(error_u,2),
        "error_p": np.linalg.norm(error_p,2),
    }

    #  Penalization Method

    t = time.time()

    eps_inv = 1/epsilon

    M = gf.Spmat('mult',Bt,B)
    M.scale(eps_inv)
    M.add(range(I_u[1]),range(I_u[1]),A)

    S = gf.linsolve('superlu',M,F)
    U_eps = np.transpose(S[0])
    P_eps = eps_inv*B*U_eps

    elapsed = time.time() - t

    error_u = u-U_eps
    error_p = p-P_eps
    #print('')
    #print('=== Penalization Method (epsilon =',epsilon,') ===')
    #print('L-2 errors w.r.t GetFEM solution (Penalization): ',np.linalg.norm(error_u,2),' / ',np.linalg.norm(error_p,2))
    #print('L-2 norm of the divergence (Penalization): ',np.linalg.norm(B*U_eps,2))
    #print('Estimation of the condition number (A+eps_inv*BtB): ', S[1])
    #print('Estimation of the memory usage (A+eps_inv*BtB): ', 8*8*3*M.nnz()*0.000008, 'Mb')
    #print('CPU Time: ', elapsed, ' s')

    tracking_penalization = {
        "method": "penalization",
        "duration": elapsed,
        "condition": S[1],
        "memory_usage": 8*8*3*M.nnz()*0.000008,
        "divergence_vector_field": np.linalg.norm(B*U_eps,2),
        "error_u": np.linalg.norm(error_u,2),
        "error_p": np.linalg.norm(error_p,2),
    }

    # Duality Method

    AinvF = np.transpose(gf.linsolve('mumps',A,F))
    BAinvF = B*AinvF

    AinvBt = gf.Spmat('empty',I_u[1],I_p[1])
    for i in range(I_p[1]):
        AinvBt.assign(range(I_u[1]),i,gf.linsolve('mumps',A,Bt[range(I_u[1]),i]))

    Uzawa = gf.Spmat('mult',B,AinvBt)

    t = time.time()

    SP = gf.linsolve('mumps',Uzawa,BAinvF)
    P_dual = np.transpose(SP[0])
    SU = gf.linsolve('superlu',A,F-Bt*P)
    U_dual = np.transpose(SU[0])

    elapsed = time.time() - t

    error_u = u-U_dual
    error_p = p-P_dual
    #print('')
    #print('=== Duality Method ===')
    #print('L-2 errors w.r.t GetFEM solution (Duality): ',np.linalg.norm(error_u,2),' / ',np.linalg.norm(error_p,2))
    #print('L-2 norm of the divergence (Duality): ',np.linalg.norm(B*U_dual,2))
    #print('Estimation of the condition numbers (A): ', SU[1])
    #print('Estimation of the memory usage (A): ', 8*8*3*A.nnz()*0.000008, 'Mb')
    #print('Estimation of the memory usage (Uzawa): ', 8*8*3*Uzawa.nnz()*0.000008, 'Mb')
    #print('CPU Time (not Uzawa computation): ', elapsed, ' s')

    tracking_duality = {
        "method": "duality",
        "duration": elapsed,
        "condition": SU[1],
        "memory_usage": 8*8*3*(A.nnz() + Uzawa.nnz())*0.000008,
        "divergence_vector_field": np.linalg.norm(B*U_dual,2),
        "error_u": np.linalg.norm(error_u,2),
        "error_p": np.linalg.norm(error_p,2),
    }

    #  Duality-Penalization Method

    Agamma = gf.Spmat('mult',Bt,B)
    Agamma.scale(Gamma)
    Agamma.add(range(I_u[1]),range(I_u[1]),A)

    AinvF = np.transpose(gf.linsolve('mumps',Agamma,F))
    BAinvF = B*AinvF

    AinvBt = gf.Spmat('empty',I_u[1],I_p[1])
    for i in range(I_p[1]):
        AinvBt.assign(range(I_u[1]),i,gf.linsolve('mumps',Agamma,Bt[range(I_u[1]),i]))

    Uzawa = gf.Spmat('mult',B,AinvBt)

    t = time.time()

    SP = gf.linsolve('mumps',Uzawa,BAinvF)
    P_dual = np.transpose(SP[0])
    SU = gf.linsolve('superlu',Agamma,F-Bt*P)
    U_dual = np.transpose(SU[0])

    elapsed = time.time() - t

    error_u = u-U_dual
    error_p = p-P_dual
    #print('')
    #print('=== Duality-Penalization Method (gamma =',Gamma,') ===')
    #print('L-2 errors w.r.t GetFEM solution (Duality-Penalization): ',np.linalg.norm(error_u,2),' / ',np.linalg.norm(error_p,2))
    #print('L-2 norm of the divergence (Duality-Penalization): ',np.linalg.norm(B*U_dual,2))
    #print('Estimation of the condition numbers (A+Gamma*BtB): ', SU[1])
    #print('Estimation of the memory usage (A+Gamma*BtB): ', 8*8*3*Agamma.nnz()*0.000008, 'Mb')
    #print('Estimation of the memory usage (Uzawa): ', 8*8*3*Uzawa.nnz()*0.000008, 'Mb')
    #print('CPU Time (not Uzawa computation): ', elapsed, ' s')

    tracking_duality_penalization = {
        "method": "duality penalization",
        "duration": elapsed,
        "condition": SU[1],
        "memory_usage": 8*8*3*(Agamma.nnz()+Uzawa.nnz())*0.000008,
        "divergence_vector_field": np.linalg.norm(B*U_dual,2),
        "error_u": np.linalg.norm(error_u,2),
        "error_p": np.linalg.norm(error_p,2), 
    }

    full_dict = {
        "method": [],
        "duration": [],
        "condition": [],
        "memory_usage": [],
        "divergence_vector_field": [],
        "error_u": [],
        "error_p": [],
        "h": [h]*4,
        "epsilon": [epsilon]*4,
        "gamma": [Gamma]*4
    }

    for dico in [tracking_global,
                 tracking_duality,
                 tracking_penalization,
                 tracking_duality_penalization]:
        for key, value in dico.items():
            full_dict[key].append(value)
    return full_dict