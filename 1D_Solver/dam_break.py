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
import shallow_roe_monthe
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
        #T0=50
        mu_r=1000
    return rho,b,cp,kappa,lambd,Tc,T0,mu_r

def qinit(state,IC='dam_break'):
    q=state.q
    x0=0.
    xc = state.grid.x.centers
    rho,b,cp,kappa,lambd,Tc,T0,mu_r=get_lava_parameters("Etna")
    T_env=300

    if IC=='dam-break':
        bath=np.zeros(np.size(xc))
        state.aux[0,:]= bath
        zl = 5.
        ul = 0.
        zr = 1.e-4
        ur = 0.
        hl= (zl-bath)*(xc<=x0)
        hr= (zr-bath)*(xc>x0)
        q[depth,:] = hl+hr
        q[momentum_U,:] = hl*ul + hr*ur
        q[momentum_T,:] = hl*T0 + hr*T0
    elif IC=='stationary-flat':
        bath=np.zeros(np.size(xc))
        state.aux[0,:]= bath
        q[depth,:] = 1- state.aux[0, :]
        q[momentum_U,:] = 0.
        q[momentum_T,:] = T0*state.q[depth,:]
    elif IC=='stationary-bump':
        bath=np.zeros(np.size(xc))
        state.aux[0,:]= state.aux[0, :] = 0.8 * np.exp(-xc**2 / 0.2**2) - 1.0
        #temp=T0*(xc>=-1)*(xc<=1)+T_env*(xc<-1)+T_env*(xc>1)
        q[depth,:] = 1- state.aux[0, :]
        q[momentum_U,:] = 0.
        q[momentum_T,:] = T0*state.q[depth,:]
    elif IC=='two-bumps':
        state.aux[0, :] = 0.8 * np.exp(-xc**2 / 0.2**2) - 1.0
        q[depth, :] = 0.1 * np.exp(-(xc + 0.4)**2 / 0.2**2) - state.aux[0, :]
        q[momentum_U, :] = 0.0
        q[momentum_T, :] = T0*q[depth, :]
    elif IC=='slope':
        bath=0.4*xc+1
        state.aux[0,:]= bath
        q[depth,:] = 5- state.aux[0, :]
        q[momentum_U,:] = 0.
        q[momentum_T,:] = T0*q[depth,:]
    elif IC=='dry-slope':
        bath=0.2*xc+1
        state.aux[0,:]=bath
        q[depth,:] = 1*(xc>=x0)+(1.e-4)*(xc<x0)
        q[momentum_U,:] = 0.
        q[momentum_T,:] = T0*q[depth,:]
    elif IC=='double-dam-break':
        bath=np.zeros(np.size(xc))
        state.aux[0,:] = bath
        q[depth,:] = 5*(xc>-0.5)*(xc<0.5)+(1.e-4)*(xc<=-0.5)*(xc>=0.5)
        q[momentum_U,:] = 0.
        q[momentum_T,:] = T0*q[depth,:]



#    elif IC=='2-shock':
#        hl = 1.
#        ul = 1.
#        hr = 1.
#        ur = -1.
#        state.q[depth,:] = hl * (xc <= x0) + hr * (xc > x0)
#        state.q[momentum,:] = hl*ul * (xc <= x0) + hr*ur * (xc > x0)
#    elif IC=='perturbation':
#        eps=0.1
#        state.q[depth,:] = 1.0 + eps*np.exp(-(xc-x0)**2/0.5)
#        state.q[momentum,:] = 0.

#Defining source term

def source_step_Lava_Flow(solver,state,dt):
    """
    Source term integration, this handles bathymetry and the temperature
    dependence of viscosity
    This is a Clawpack-style source term routine, which approximates
    the integral of the source terms over a step.
    """
    #Obtaining data from state
    dryt=state.problem_data['dry_tolerance'] 
    q = state.q
    bath=state.aux
    xc = state.grid.x.centers
    
    #Obtaining parameters characteristic to the lava
    rho,b,cp,kappa,lambd,Tc,T0,mu_r=get_lava_parameters("Etna")
    nu_r=mu_r/rho
    T_env=300
    
    #Computing slope for bathymetry
    dx=xc[3]-xc[2]
    dHdx=np.zeros(np.size(bath[0]))
    dHdx[1:]=(1/dx)*(bath[0,1:]-bath[0,0:-1])
    dHdx[0]=dHdx[1]
 
    h = q[0,:]
    u   = q[1,:]/h
    T   = q[2,:]/h

    #Computing more parameters (these are for Etna lava flows)
    H=(3*1.e-6)/h
    E=1.5*1.e-15
    W=2*1.e-6
    K=(4*1.e-3)/h
    
    #Integrating one step in time just of the cell is not dry
    q[1,:] = np.where(h<=dryt,  q[1,:],  (h**2)*q[1,:]/(h**2+3*nu_r*dt*np.exp(-b*(T-T0)))-grav*dt*h*dHdx)
    q[2,:] = np.where(h<=dryt, q[2,:], q[2,:]+dt*(-E*(T**4-T_env**4)-W*(T-T_env)-H*(T-Tc)+K*(u**2)*np.exp(-b*(T-T0))))

def setup(use_petsc=False,kernel_language='Fortran',outdir='./_output',solver_type='classic',
          riemann_solver='monthe', disable_output=False):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language == 'Python':
        if riemann_solver.lower() == 'roe':
            raise Exception('Python Roe solver not implemented.')
        elif riemann_solver.lower() == 'hlle':
            raise Exception('Python HLLE solver not implemented.')
    elif kernel_language == 'Fortran':
        if riemann_solver.lower() == 'monthe':
            rs = shallow_roe_monthe
 
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(rs)
        solver.order=1
        #solver.limiters = pyclaw.limiters.tvd.vanleer
        #Adding source term to conservation laws
        solver.step_source=source_step_Lava_Flow
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)

    solver.kernel_language = kernel_language

    solver.bc_lower[0] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.extrap
    #Auxiliary vector will contain bathymetry
    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap

    solver.num_eqn=num_eqn
    solver.num_waves=num_waves
   
    xlower = lowerx
    xupper = upperx
    mx = 500
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,num_eqn,num_aux=1)

    # Gravitational constant
    state.problem_data['grav'] = 9.8
    state.problem_data['dry_tolerance'] = 1e-3
    state.problem_data['sea_level'] = 0.0
    

    IC='dry-slope'
    qinit(state,IC=IC)


    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.tfinal = 30.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
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
    plotaxes.title = 'Water height'
    plotaxes.axescmd = 'subplot(311)'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta
    plotitem.plot_var2= bathy
    plotitem.color = rgb_converter((67,183,219))
    
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
