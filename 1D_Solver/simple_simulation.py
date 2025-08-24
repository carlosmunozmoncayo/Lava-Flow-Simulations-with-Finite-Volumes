#!/usr/bin/env python
# encoding: utf-8

r"""
Shallow water flow
==================

Solve the one-dimensional shallow water equations:

.. math::
    h_t + (hu)_x & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x & = 0.

Here h is the depth, u is the velocity, and g is the gravitational constant.
The default initial condition used here models a dam break.
"""

from __future__ import absolute_import
import numpy as np
from clawpack import riemann
import geoclaw_swe_rs
import newton_solver

#from clawpack.riemann.shallow_roe_with_efix_1D_constants import depth, momentum, num_eqn
depth=0
momentum_U=1
momentum_T=2
num_eqn=3
num_waves=3
grav=9.8

lowerx=-5.0
upperx=5.0


def get_lava_parameters(lava_flow):
    if lava_flow=="Etna":
        rho=2500
        b=0.02
        cp=1200
        kappa=2.0
        lambd=70
        Tc=1253
        T0=1353
        mu_r=1000
    return rho,b,cp,kappa,lambd,Tc,T0,mu_r

def qinit(state,IC='dam_break',Tvent=1200):
    q=state.q
    dryt=state.problem_data['dry_tolerance'] 
    x0=0.
    xc = state.grid.x.centers
    rho,b,cp,kappa,lambd,Tc,T0,mu_r=get_lava_parameters("Etna")
    T_env = 300
    T0 = Tvent #Cold lava

    if IC=='dam-break':
        bath=np.zeros(np.size(xc))
        state.aux[0,:]= bath
        ul = 0.
        ur = 0.
        hl= 2.
        hr= 1.e-4
        q[depth,:] = hl*(xc<=x0)+hr*(xc>x0)
        q[momentum_U,:] = hl*ul + hr*ur
        q[momentum_T,:] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])
    elif IC=='dam-break-sharp-obstacle':
        bath=2*(xc>1)*(xc<2)
        state.aux[0,:]=bath 
        ul = 0.
        ur = 0.
        hl= 2
        hr= 1.e-4
        q[depth,:] = hl*(xc<=x0)+hr*(xc>x0)
        q[momentum_U,:] = 0.
        q[momentum_T,:] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])
    elif IC=='dam-break-smooth-obstacle':
        bath=2 * np.exp(-(xc-2)**2 *10) 
        state.aux[0,:]=bath
        zl = 5.0
        ul = 0.
        zr = 1.e-4
        ur = 0.
        hl= (zl-bath)*(xc<=x0)
        hr= (zr-bath)*(xc>x0)
        q[depth,:] = 2*(xc<=x0)+(1.e-4)*(xc>x0)
        q[momentum_U,:] = 0 
        q[momentum_T,:] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])

    elif IC=='stationary-flat':
        bath=np.zeros(np.size(xc))
        state.aux[0,:]= bath
        q[depth,:] = 1
        q[momentum_U,:] = 0.
        q[momentum_T,:] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])
    elif IC=='stationary-bump':
        bath=np.zeros(np.size(xc))
        state.aux[0,:]= 0.8 * np.exp(-xc**2 / 0.2**2)
        #temp=T0*(xc>=-1)*(xc<=1)+T_env*(xc<-1)+T_env*(xc>1)
        q[depth,:] = 1- state.aux[0, :]
        q[momentum_U,:] = 0.
        q[momentum_T,:] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])
    elif IC=='stationary-circle':
        bath=xc**2/4
        state.aux[0,:]=bath
        q[depth,:] = np.where(2-bath>=0, 2-bath, 1.e-4)
        q[momentum_U,:] = 0.
        q[momentum_T,:] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])

    elif IC=='two-bumps':
        state.aux[0, :] = 0.8 * np.exp(-xc**2 / 0.2**2) - 1.0
        q[depth, :] = 0.1 * np.exp(-(xc + 0.4)**2 / 0.2**2) - state.aux[0, :]
        q[momentum_U, :] = 0.0
        q[momentum_T, :] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])
    elif IC=='slope':
        bath=0.4*xc+1
        state.aux[0,:]= bath
        q[depth,:] = 5- state.aux[0, :]
        q[momentum_U,:] = 0.
        q[momentum_T,:] =  np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])
    elif IC=='dry-slope':
        bath=0.2*xc+1
        state.aux[0,:]=bath
        q[depth,:] = 1*(xc>=x0)+(1.e-4)*(xc<x0)
        q[momentum_U,:] = 0.
        q[momentum_T,:] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])
    elif IC=='dry-slope-sharp-obstacle':
        bath=0.5*xc+2*(xc>-2)*(xc<-1)
        state.aux[0,:]=bath
        q[depth,:] = 3*(xc>=x0)+(1.e-4)*(xc<x0)
        q[momentum_U,:] = 0.
        q[momentum_T,:] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])
    elif IC=='dry-slope-smooth-obstacle':
        bath=0.5*xc+2 * np.exp(-(xc+2)**2 *10)
        state.aux[0,:]=bath
        q[depth,:] = 3*(xc>=x0)+(1.e-4)*(xc<x0)
        q[momentum_U,:] = 0.
        q[momentum_T,:] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])
    elif IC=='double-dam-break':
        bath=np.zeros(np.size(xc))
        state.aux[0,:] = bath
        q[depth,:] = 5*(xc>-0.5)*(xc<0.5)+(1.e-4)*(xc<=-0.5)*(xc>=0.5)
        q[momentum_U,:] = 0.
        q[momentum_T,:] = np.where(q[depth,:]<=dryt, T_env*q[depth,:],T0*q[depth,:])
    elif IC=='volcano-column':
        slope=1
        bath=(slope*xc+2)*(xc<=0)*(xc>-2)-(slope*xc-2)*(xc>0)*(xc<2)
        state.aux[0,:] = bath
        q[depth,:] = (4+slope-bath)*(xc>-0.5)*(xc<0.5)+(1.e-4)*(xc<=-0.5)*(xc>=0.5)
        q[momentum_U,:] = 0.
        q[momentum_T,:] = T0*q[depth,:]

def source_step_Lava_Flow(solver,state,dt):
    """
    Source term integration, this handles bathymetry and the temperature
    dependence of viscosity
    This is a Clawpack-style source term routine, which approximates
    the integral of the source terms over a step.
    """
    #Obtaining data from state
    dryt= state.problem_data['dry_tolerance'] 
    friction_coeff= state.problem_data['friction_coeff']
    friction_depth = state.problem_data['friction_depth']
    g= state.problem_data['grav']
    q = state.q
    bath=state.aux
    xc = state.grid.x.centers

    h   = q[0,:]
    u   =  np.where(h<=dryt, 0, q[1,:]/h)

    #Semi-implicit integration of friction source term
    if friction_coeff>0:
        logical_var = np.logical_and(h>dryt, h<=friction_depth)
        gamma = np.where(logical_var, 0, np.abs(q[1,:])*(g*friction_coeff**2)/(h**(7.0/3.0)))
        q[1,:] = np.where(logical_var, 0, q[1,:]/(1.0 + dt*gamma))
    
    #Obtaining parameters characteristic to the lava
    rho,b,cp,kappa,lambd,Tc,T0,mu_r=get_lava_parameters("Etna")
    nu_r=mu_r/rho
    T_env=300
    T   = np.where(h<=dryt, T_env, q[2,:]/h)

    #Computing more parameters (these are for Etna lava flows)
    H=np.where(h<=dryt, 0, (3*1.e-6)/h)
    E=1.5*1.e-15
    W=2*1.e-6
    K=np.where(h<=dryt, 0, (4*1.e-3)/h)
    
    #Integrating one step in time just if the cell is not dry
    q[1,:] = np.where(h<=dryt,  q[1,:],  (h**2)*q[1,:]/(h**2+3*nu_r*dt*np.exp(-b*(T-T0))))
    q[2,:] = np.where(h<=dryt,  q[2,:], q[2,:]+dt*(-E*(T**4-T_env**4)-W*(T-T_env)-H*(T-Tc)+K*(u**2)*np.exp(-b*(T-T0))))

def source_step_Lava_Flow_implicit(solver,state,dt):
    #Taking a backward Euler step with a Newton solver for the source term
    q = state.q
    drytol = state.problem_data['dry_tolerance']
    g = state.problem_data['grav']

    qp = q.copy() #Primitive variables
    qp[1,:] = np.where(q[0,:]<=drytol, 0, q[1,:]/q[0,:]) #Velocity
    qp[2,:] = np.where(q[0,:]<=drytol, 300, q[2,:]/q[0,:]) #Temperature

    qnp1,status,iters=newton_solver.newton_solver_mod.newton_batch(qn=qp,dt=dt,tol_step=1e-14,tol_f=1.e-14,maxit=50)
    print(np.max(np.abs(qnp1-qp)))
    print(np.max(iters))
    state.q[1,:] = np.where(q[0,:]<=drytol, q[1,:], qnp1[1,:]*q[0,:])
    state.q[2,:] = np.where(q[0,:]<=drytol, q[2,:], qnp1[2,:]*q[0,:])

def clean_dry_cells_b4step(solver,state):
    """
    This function cleans the dry cells before each time step
    """
    dryt=state.problem_data['dry_tolerance'] 
    q = state.q
    h = q[0,:]
    #Cleaning the dry cells
    # state.q[0,:] = np.where(h<=dryt, 0, q[0,:])
    state.q[1,:] = np.where(h<=dryt, 0., q[1,:])
    state.q[2,:] = np.where(h<=dryt, 0., q[2,:]) #Environmental temperature

def setup(use_petsc=False,kernel_language='Fortran',outdir='./_output',solver_type='classic',
        riemann_solver='monthe', disable_output=False, IC='dry-slope-smooth-obstacle',
        nframes=10, tfinal=5.0,mx=500,Tvent=1200):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language == 'Python':
            raise Exception('Python Riemann solver not implemented.')
    elif kernel_language == 'Fortran':
            rs = geoclaw_swe_rs
 
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(rs)
        solver.order=1
        #solver.limiters = pyclaw.limiters.tvd.vanleer
        #Adding source term to conservation laws
        # solver.step_source=source_step_Lava_Flow
        solver.step_source=source_step_Lava_Flow_implicit
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)

    solver.kernel_language = kernel_language
    solver.fwave = True

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.wall
    #Auxiliary vector will contain bathymetry
    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap

    solver.num_eqn=num_eqn
    solver.num_waves=num_waves

    solver.cfl_max=0.95
    solver.cfl_desired=0.85

    # solver.before_step = clean_dry_cells_b4step
   
    xlower = lowerx
    xupper = upperx
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,num_eqn,num_aux=1)
    
    ###############################################################
    ############# Problem parameters ##############################
    ###############################################################
    # Gravitational constant, water, and ground parameters
    state.problem_data['grav'] = 9.81
    state.problem_data['dry_tolerance'] = 1e-3
    state.problem_data['sea_level'] = 0.0
    state.problem_data['friction_coeff'] = 0.025 #Manning coefficient for rocky bed (?)
    state.problem_data['friction_depth'] = 0.1 #Depth for friction to act (m)

    #Lava flow parameters
    rho,b,cp,kappa,lambd,Tc,T0,mu_r=get_lava_parameters("Etna")
    nur=mu_r/rho
    Tenv=300
    params_in = [nur, b, Tc, Tenv, T0,
                  3*1.e-6, 4*1.e-3, 2*1.e-6, 1.5*1.e-15,
                  state.problem_data['dry_tolerance'], state.problem_data['grav']]
    newton_solver.newton_solver_mod.set_params(params_in)
    newton_solver.newton_solver_mod.print_params()
    ###############################################################
    ###############################################################

    

    qinit(state,IC=IC,Tvent=Tvent)


    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.tfinal = tfinal
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.num_output_times = nframes
    claw.setplot = setplot
    claw.write_aux_init = True

    return claw


#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data



    # Plot variables
    def bathy(current_data):
        return current_data.aux[0, :]

    def eta(current_data):
        return current_data.q[0, :] + bathy(current_data)

    def velocity(current_data):
        return current_data.q[1, :] / current_data.q[0, :]

    def temperature(current_data):
        return current_data.q[2, :] / current_data.q[0, :]

    rgb_converter = lambda triple: [float(rgb) / 255.0 for rgb in triple]

    # Figure for depth
    plotfigure = plotdata.new_plotfigure(name='Lava height', figno=0)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [lowerx,upperx]
    plotaxes.ylimits = [-3,5]#[-0.1,0.5]#
    plotaxes.title = 'Water height'
    plotaxes.axescmd = 'subplot(311)'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta
    plotitem.plot_var2= bathy
    plotitem.color = rgb_converter((255, 69, 0))
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = bathy
    plotitem.color = 'k'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta
    plotitem.color = 'k'

    # Figure for momentum[1]
    #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = [lowerx,upperx]
    plotaxes.title = 'Velocity'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = velocity#momentum_U
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':1}
    
    #Figure for momentum[2]
    #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = [lowerx,upperx]
    plotaxes.title = 'Temperature'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = temperature#momentum_T
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':1}
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
