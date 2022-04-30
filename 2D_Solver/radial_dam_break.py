#!/usr/bin/env python
# encoding: utf-8
r"""
2D Lava flow: radial dam break
==================================

Solve the 2D lava flow equations proposed by Costa and Macedonio (2005):

.. math::
    h_t + (hu)_x + (hv)_y = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x + (huv)_y = -g h H_x-\gamma u \\
    (hv)_t + (huv)_x + (hv^2 + \frac{1}{2}gh^2)_y = -g h H_y-\gamma v \\
    (hT)_t + (huT)_x + (hvT)_y = -\mathcal{E}\left(T^{4}-T_{e n v}^{4}\right)-\mathcal{W}\left(T-T_{e n v}\right)\\
    -\mathcal{H}\left(T-T_{c}\right)+\mathcal{K}& \left(U^{2}+V^{2}\right) \exp \left[-b\left(T-T_{r}\right)\right]


The initial condition is a circular area with high depth surrounded by lower-depth lava.
The top and right boundary conditions reflect, while the bottom and left boundaries
are outflow.
"""

from __future__ import absolute_import
import numpy as np
import shallow_roe_monthe
from clawpack import riemann
import sys

#We add a four equation for temperature
num_eqn=4

depth=0
x_momentum=1
y_momentum=2
T_momentum=3

def qinit(state,h_in=4.,h_out=1.,dam_radius=0.5):
    x0=0.
    y0=0.
    X, Y = state.p_centers
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    rho,b,cp,kappa,lambd,Tc,T0,mu_r=get_lava_parameters("Etna")

    state.q[depth     ,:,:] = h_in*(r<=dam_radius) + h_out*(r>dam_radius)
    state.q[x_momentum,:,:] = 0.
    state.q[y_momentum,:,:] = 0.
    state.q[T_momentum,:,:] = T0*state.q[depth     ,:,:]   #T0*h_in*(r<=dam_radius) + Tc*h_out*(r>dam_radius)

#The model requires several parameters that depend on the lava we are trying to model
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
    for i in range(11):
        q = state.q

        h = q[0,:,:]
        u   = q[1,:,:]/h
        v   = q[2,:,:]/h
        T  = q[3,:,:]/h
        
        #Computing more parameters (these are for Etna lava flows)

        H=(3*1.e-6)/h
        E=1.5*1.e-15
        W=2*1.e-6
        K=(4*1.e-3)/h

        #q[0,:,:] = q[0,:,:] - dt/rad * qstar[2,:,:]
        q[1,:,:] = (h**2)*q[1,:,:]/(h**2+3*nu_r*ddt*np.exp(-b*(T-T0)))
        q[2,:,:] = (h**2)*q[2,:,:]/(h**2+3*nu_r*ddt*np.exp(-b*(T-T0)))
        q[3,:,:] = q[3,:,:]+dt*(-E*(T**4-T_env**4)-W*(T-T_env)-H*(T-Tc)+K*(u**2+v**2)*np.e**(-b*(T-T0)))



def setup(kernel_language='Fortran', use_petsc=False, outdir='./_output',
          solver_type='classic', riemann_solver='monthe',disable_output=False):
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if riemann_solver.lower() == 'roe':
        rs = riemann.shallow_roe_with_efix_2D
    elif riemann_solver.lower() == 'hlle':
        rs = riemann.shallow_hlle_2D
    elif riemann_solver.lower() == 'monthe':
        rs = shallow_roe_monthe


    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(rs)
        solver.limiters = pyclaw.limiters.tvd.MC
        solver.dimensional_split=True
        #If the line below is commented, just the homogeneous hyperbolic system is being solved
        #solver.step_source=source_step_Lava_Flow
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(rs)

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.wall
    solver.bc_lower[1] = pyclaw.BC.extrap
    solver.bc_upper[1] = pyclaw.BC.wall

    # Domain:
    xlower = -2.5
    xupper = 2.5
    mx = 150
    ylower = -2.5
    yupper = 2.5
    my = 150
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])

    state = pyclaw.State(domain,num_eqn)


    # Gravitational constant
    state.problem_data['grav'] = 9.8

    qinit(state)

    claw = pyclaw.Controller()
    claw.tfinal = 0.1
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    solver.num_eqn=4
    solver.num_waves=4
    #print("num_waves: ", solver.num_waves)
    if disable_output:
        claw.output_format = None
    claw.outdir = outdir
    claw.num_output_times = 10
    claw.setplot = setplot
    claw.keep_copy = True

    return claw

#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for depth
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2.5, 2.5]
    plotaxes.ylimits = [-2.5, 2.5]
    plotaxes.title = 'Water height'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    #plotitem.pcolor_cmin = 0.5
    #plotitem.pcolor_cmax = 1.5
    plotitem.add_colorbar = True
    
    # Scatter plot of depth
    plotfigure = plotdata.new_plotfigure(name='Scatter plot of h', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0., 2.5]
    plotaxes.ylimits = [0., 2.1]
    plotaxes.title = 'Scatter plot of h'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = depth
    def q_vs_radius(current_data):
        from numpy import sqrt
        x = current_data.x
        y = current_data.y
        r = sqrt(x**2 + y**2)
        q = current_data.q[depth,:,:]
        return r,q
    plotitem.map_2d_to_1d = q_vs_radius
    plotitem.plotstyle = 'o'


    # Figure for x-momentum
    plotfigure = plotdata.new_plotfigure(name='Momentum in x direction', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2.5, 2.5]
    plotaxes.ylimits = [-2.5, 2.5]
    plotaxes.title = 'Momentum in x direction'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = x_momentum
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = False       # show on plot?
    

    # Figure for y-momentum
    plotfigure = plotdata.new_plotfigure(name='Momentum in y direction', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2.5, 2.5]
    plotaxes.ylimits = [-2.5, 2.5]
    plotaxes.title = 'Momentum in y direction'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = y_momentum
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = False       # show on plot?
    


    # Figure for T-momentum
    plotfigure = plotdata.new_plotfigure(name='Momentum  Temperature', figno=4)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2.5, 2.5]
    plotaxes.ylimits = [-2.5, 2.5]
    plotaxes.title = 'Momentum Temperature'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = T_momentum
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
