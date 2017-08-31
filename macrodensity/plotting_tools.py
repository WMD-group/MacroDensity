import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def energy_band_alignment_diagram(energies,materials,limit=8.,width=1.,cols=['#74356C','#efce19'],textsize=22,arrowhead=0.7, outfile='BandAlignment', references=[]):
    '''
    Function for plotting the classic energy band alignment diagram
    Args:
	energies: a list of EAs and IPs for all materials in the plot [[ea1,ip1],[ea2,ip2], ....]
	materials: a list of labels for the materials
	limit: the deepest value you want on the plot (eg 8 means that you can plot IPs up to 8, in practive you want this higher than the greatest IP, usually by about 2 eV). Default = 8.
	width: The width of the bars, nearly always 1. Default = 1.
	cols: The list of colours that you want to rotate through. Default = [plum, buttercup]
	textsize: size of the font for the figure. Default = 22.
	arrowhead: arrow head length. Default = 0.7.
	outfile: name of the ouput. Default = 'BandAlignment'
        references: any reference levels you want to add to the plot. [['Name of reference',value_of_reference], ...]. Note that value_of_reference is a positive value on the same scale as IP/EA. Default = [].
    Returns:
        Nothing, but draws an eps plot.
    '''

    fig, ax1 = plt.subplots(1, 1, sharex=True)
    fig.set_size_inches(len(energies) * 3, limit * 0.75)
    mpl.rcParams['xtick.labelsize'] = textsize
    mpl.rcParams['ytick.labelsize'] = textsize
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['ytick.major.width'] = 3
    mpl.rcParams['ytick.major.size'] = 7
    mpl.rcParams['ytick.minor.size'] = 4
    mpl.rcParams['axes.linewidth'] = 3

    ax1.set_color_cycle(cols)
    ax2 = ax1.twinx()
    ind = np.arange(len(energies))

## Bars for the IP and background colour
    for i in ind:
        ax1.bar(i,-limit, width, edgecolor='black')
        ax1.bar(i,-energies[i][1], width, color='w', edgecolor='black')
    
## Reset the colours back to the start and plot the EA
    ax1.set_color_cycle(cols)
    for i in ind:
        ax1.bar(i,-energies[i][0], width, edgecolor='black',alpha=0.8)

## Set the limits of the axes
    ax1.set_ylim(-limit,0)
    ax2.set_ylim(-limit,0)
    ax1.set_xlim(-0.5,len(energies)-0.5)

## Set the names
    ax1.set_xticks(ind)
    ax1.set_xticklabels(materials,size=textsize)
    ran = [ str(k) for k in np.arange(0,limit+2,2)]
    ax1.set_yticklabels(ran[::-1],size=textsize)
    ran = [ '' for k in np.arange(0,limit+2,2)]
    ran[0] = 'Vacuum Level'
    ax2.set_yticklabels(ran[::-1],size=textsize)
    ax1.set_ylabel('Energy (eV)', size=textsize)


    os1 = 0.15   # Offset of the text 'IP' in the plot
    os2 = 0.2    # Offset of the text 'EA' in the plot

    for i, en in enumerate(energies):
        ax1.arrow(i,-en[0],0, en[0]-arrowhead, width=0.005, head_length=arrowhead,head_width=0.07, fc='black',ec='None')
        ax1.arrow(i,0,0, -en[1]+arrowhead, width=0.005, head_length=arrowhead,head_width=0.07, fc='black',ec='None')
        ax1.arrow(i,0,0, -en[0]+arrowhead, width=0.005, head_length=arrowhead,head_width=0.07, fc='black',ec='None')
        loc_ip = -(en[0] + en[1]) / 2
        ax1.text(i+os1*0.1,loc_ip,en[1],fontsize=textsize)
        ax1.text(i-os1,loc_ip,'IP',fontsize=textsize)
    
        loc_ea = -en[0] / 2
        ax1.text(i+os1*0.1,loc_ea,en[0],fontsize=textsize)
        ax1.text(i-os2,loc_ea,'EA',fontsize=textsize)
    
        ax1.minorticks_on()
        ax1.tick_params(axis='x',which='minor',bottom='off') # Don't show minor ticks on x-axis
        ax2.minorticks_on()

    for ref in references:
        ax1.hlines(-ref[1],-0.5,len(energies)-0.5,linestyles='--',colors='r')
        ax1.text(len(energies)-0.45,-ref[1]-0.1,ref[0],fontsize=textsize,color='r')

    fig.savefig('%s.eps'%outfile,bbox_inches='tight')
    print "Figure saved as %s.eps"%(outfile)


