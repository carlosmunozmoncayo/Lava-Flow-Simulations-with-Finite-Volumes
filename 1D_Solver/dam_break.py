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
    x0=0.
    xc = state.grid.x.centers
    rho,b,cp,kappa,lambd,Tc,T0,mu_r=get_lava_parameters("Etna")
    if IC=='dam-break':
        hl = 3.
        ul = 0.
        hr = 1.
        ur = 0.
        state.q[depth,:] = hl * (xc <= x0) + hr * (xc > x0)
        state.q[momentum_U,:] = hl*ul * (xc <= x0) + hr*ur * (xc > x0)
        state.q[momentum_T,:] = hl*T0 * (xc <= x0) + hr*T0 * (xc > x0)
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
    #Obtaining sparameters ccharacteristic to the lava
    rho,b,cp,kappa,lambd,Tc,T0,mu_r=get_lava_parameters("Etna")
    nu_r=mu_r/rho
    T_env=300

    #Doing several time integrations of source term to avoid stiffness
    ddt=dt/10
    if True:
        q = state.q

        h = q[0,:]
        u   = q[1,:]/h
        T   = q[2,:]/h
        #T  = q[3,:,:]/h

        #Computing more parameters (these are for Etna lava flows)

        H=(3*1.e-6)/h
        E=1.5*1.e-15
        W=2*1.e-6
        K=(4*1.e-3)/h

        #q[0,:,:] = q[0,:,:] - dt/rad * qstar[2,:,:]
        q[1,:] = (h**2)*q[1,:]/(h**2+3*nu_r*dt*np.exp(-b*(T-T0)))
        #q[2,:] = (h**2)*q[2,:]/(h**2+3*nu_r*ddt*np.exp(-b*(T-T0)))
        q[2,:] = q[2,:]+dt*(-E*(T**4-T_env**4)-W*(T-T_env)-H*(T-Tc)+K*(u**2)*np.exp(-b*(T-T0)))
    else:
        q = state.q
        q[2,:]=4000

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
        #solver.step_source=source_step_Lava_Flow
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)

    solver.kernel_language = kernel_language

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap

    solver.num_eqn=num_eqn
    solver.num_waves=num_waves

    xlower = -5.0
    xupper = 5.0
    mx = 500
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,num_eqn)

    # Gravitational constant
    state.problem_data['grav'] = 1.0
    state.problem_data['dry_tolerance'] = 1e-3
    state.problem_data['sea_level'] = 0.0
    

    IC='dam-break'
    qinit(state,IC=IC)


    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.tfinal = 2.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot

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

    # Figure for depth
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Water height'
    plotaxes.axescmd = 'subplot(311)'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = depth
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}

    # Figure for momentum[1]
    #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Momentum'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = momentum_U
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}
    
    #Figure for momentum[2]
    #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Momentum_T'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = momentum_T
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
