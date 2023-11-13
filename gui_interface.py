# -*- coding: utf-8 -*-
"""
gui_interface.py
Interface for USC workshop
Created on Mon Sep 18 20:45:50 2023
@author:Carlos Alberto Duran Villalobos-The University of Manchester
"""
import ipywidgets as widgets
from IPython.display import display, clear_output, Javascript, HTML
from IPython.core.display import display, HTML
from cells_sim_2 import tcellmodel, control_tcell
import numpy as np
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn import preprocessing
import ctypes


def sim(tcellmodel, hours, us, variability):
    #plt.style.available
    plt.style.use('seaborn-darkgrid')#seaborn-v0_8') #use plot styles
    [sx,sy,su] = control_tcell(tcellmodel, hours, us, variability)
    h_cells = sy[96] * sx[96, 3]
    
    #plot figures
    font = 10
    fig, axs = plt.subplots(1,2, figsize=(19 / 2.54, 6/ 2.54))
    
    axs[0].plot(su,color='tab:brown', linewidth=3)
    axs[0].set_xlabel('Time [h]', fontsize = font)
    axs[0].set_ylabel('Feed rate [mL/h]', fontsize = font,color='tab:brown')
    axs[0].yaxis.set_tick_params(labelsize=font)
    axs[0].xaxis.set_tick_params(labelsize=font);
    axs[0].tick_params(axis='y', labelcolor='tab:brown')
    axs[0].set_ylim([0, 2])
    axs[1].plot(sy*sx[:,3],color='blue',linewidth=3, label='Viable cells')
    
    x=np.array([0.1,0.1,0.3,0.5,1.0,0.8,0.7,0.1,0.1,0.1])
    us=np.array([])
    for j in range(10):
        us = np.append(us, x[j]*np.ones(12)) if us.size else x[j]*np.ones(12)
    [sx,sy,su] = control_tcell(tcellmodel, hours, us, variability)
    axs[1].plot(sy*sx[:,3],color='blue',alpha=0.3 , linewidth=3, label='Target')
    
    #axs[1].plot(sy,color='blue', linewidth=3)
    axs[1].set_xlabel('Time [h]', fontsize = font)
    axs[1].set_ylabel('Viable T-cells', fontsize = font,color='blue')
    axs[1].axhline(y=100e6, color="red", linestyle="--", label='Harvest')
    axs[1].axvline(x=96, color="red", linestyle="--", label='Harvest')
    axs[1].yaxis.set_tick_params(labelsize=font)
    axs[1].xaxis.set_tick_params(labelsize=font)
    #
    axs[1].tick_params(axis='y', labelcolor='blue')
    axs[1].set_ylim([0, 1.7e8])
    handles, labels = axs[1].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    axs[1].legend(by_label.values(), by_label.keys(), fontsize = font)

    #fig.delaxes(axs[2,2])
    fig.tight_layout(pad=0.8)
    clear_output(wait = True)
    #plt.pause(0.1)
    plt.show()
    
    return h_cells


def display_selected_plot(a,b,c,d,e,f,g,h,i,j,via,inh,rg,bi):
        #Simulation with a feed
    hours=119
    #var = [95, -0.5, 0.5, -0.5] #variability: [viability, inhibitor, rg, bi]
    var = [via/100, inh/100, rg/100, bi/100]
    x=np.array([a,b,c,d,e,f,g,h,i,j])
    us=np.array([])
    for j in range(10):
        us = np.append(us, x[j]*np.ones(12)) if us.size else x[j]*np.ones(12)
    cells = sim(tcellmodel, hours, us, var)
    
    bar=widgets.FloatProgress(
                value=sum(us),
                min=0,
                max=238.0,
                description='Total feed ml:',
                bar_style='info',
                style={'bar_color': 'lightcoral'},
                orientation='horizontal'
    )
    
    
    # Calculate the value to be displayed
    value_to_display = sum(us)

    # Create widgets for the components of the text with styling for the mL label
    label_1 = widgets.HTML(value=f"<span style='font-size: 12pt; font-weight: bold'>{value_to_display:.2f}</span>")
    label_2 = widgets.HTML(value="<span style='font-size: 12pt; font-weight: bold'>mL</span>")

    # Create an HBox to arrange the widgets
    text = widgets.HBox([label_1, label_2])
    

    # Create widgets for the components of the text
    label_1 = widgets.HTML(value="<span style='font-size: 12pt'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Viable T-cells:</span>")
    label_2 = widgets.HTML(value=f"<span style='font-size: 12pt; font-weight: bold'>{cells:.4e}</span>")

    # Create an HBox to arrange the widgets
    text_cells = widgets.HBox([label_1, label_2])

    box_layout3 = widgets.Layout(display='flex',
                        flex_flow='row',
                        align_items='stretch',
                        border='solid',
                        width='100%')
    
    items = [bar, text, text_cells]
    box2 = widgets.Box(children=items, layout=box_layout3)
    user32 = ctypes.windll.user32
    dpi = user32.GetDpiForSystem()
    box2.layout.width = f'{int(20/ 2.54 * dpi)}px'
    display(box2)
    
    
    
def gluc_slider():
    user32 = ctypes.windll.user32
    dpi = user32.GetDpiForSystem()
    grid = widgets.GridspecLayout(1, 10)
    time =['0-12 h','12-24 h','24-36 h','36-48 h','48-60 h','60-72 h','72-84 h','84-96 h','96-108 h','108-120 h']
    
    for i in range(1):
        for j in range(10):
            grid[i, j] = widgets.FloatSlider(value=0.5,min=0,max=1,step=0.1,orientation='vertical',description=time[j],handle_style = {'size':20}) 
            grid[i,j].layout.width = f'{int(1.7/ 2.54 * dpi)}px'
            grid[i,j].layout.height = f'{int(4.3/ 2.54 * dpi)}px'
            grid[i,j].style.handle_color = 'lightcoral'  # Change the handle color
            grid[i,j].style.bar_color = 'lightblue'  # Change the bar color
            
    h=widgets.HTML(value='<{size}>Glucose feed rate [mL/h]:</{size}>'.format(size='h4'))
    h.layout.width = f'{int(18/ 2.54 * dpi)}px'
    h.layout.height = f'{int(0.8/ 2.54 * dpi)}px'
    
    box_layout2 = widgets.Layout(display='flex',
                        flex_flow='column',
                        align_items='stretch',
                        border='solid',
                        width='100%')
    
    grid.layout.width = f'{int(18/ 2.54 * dpi)}px'
    grid.layout.height = f'{int(4.6/ 2.54 * dpi)}px'
    items = [h, grid]
    box2 = widgets.Box(children=items, layout=box_layout2)
    box2.layout.width = f'{int(20/ 2.54 * dpi)}px'
    box2.layout.height = f'{int(5.6/ 2.54 * dpi)}px'
    box2
    
    #hide_code()
    return box2, grid

def sim_ideal(var):
    [box, grid] = gluc_slider()
    hidden = hidden_widget(var)
    out = widgets.interactive_output(display_selected_plot, {'a':grid[0,0],'b':grid[0,1],'c':grid[0,2],'d':grid[0,3],'e':grid[0,4],'f':grid[0,5],'g':grid[0,6],'h':grid[0,7],'i':grid[0,8],'j':grid[0,9],
                                                             'via':hidden[0], 'inh':hidden[1],  'rg':hidden[2], 'bi':hidden[3]})
    display(box, out)

def hidden_widget(var):
    hidden=[]
    hidden.append(widgets.IntText(value=var[0], description='Hidden Value:', style={'description_width': '0px'}))
    hidden[0].layout.display = 'none'  # Hide the widget
    hidden.append(widgets.IntText(value=var[1], description='Hidden Value:', style={'description_width': '0px'}))
    hidden[1].layout.display = 'none'  # Hide the widget
    hidden.append(widgets.IntText(value=var[2], description='Hidden Value:', style={'description_width': '0px'}))
    hidden[2].layout.display = 'none'  # Hide the widget
    hidden.append(widgets.IntText(value=var[3], description='Hidden Value:', style={'description_width': '0px'}))
    hidden[3].layout.display = 'none'  # Hide the widget
    return hidden
  
