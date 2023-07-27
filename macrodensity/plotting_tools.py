from __future__ import division, print_function

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def energy_band_alignment_diagram(energies, materials, limit=8., width=1.,
                                  cols=['#74356C','#efce19'], textsize=22,
                                  arrowhead=0.7, outfile='BandAlignment',
                                  references=[], edge=None):
   
    """
    Plot an energy band alignment diagram for a list of materials.

    Parameters:
        energies (list): A list of tuples containing the ionization potential (IP) and
                         electron affinity (EA) of each material. The format is [(IP_1, EA_1), ...].
        materials (list): A list of material names corresponding to each set of energies.
        limit (float, optional): The limit for the energy axis (in eV). Default is 8.0.
        width (float, optional): The width of the bars representing IP and EA. Default is 1.0.
        cols (list, optional): A list of colors to use for the bars. Default is ['#74356C','#efce19'].
        textsize (int, optional): The font size for the text in the plot. Default is 22.
        arrowhead (float, optional): The size of the arrowhead for the energy arrows. Default is 0.7.
        outfile (str, optional): The base name for the output file (both .eps and .png files will be saved).
                                 Default is 'BandAlignment'.
        references (list, optional): A list of reference points (as tuples) to be shown as dashed lines
                                     on the plot. Each tuple should be in the format (reference_value, label).
                                     Default is an empty list.
        edge (None or str, optional): The edge color for the bars. If None, there will be no edge color.
                                      Default is None.

    Returns:
        None: The function generates and displays the energy band alignment diagram.

    Example:
        >>> energies = [(5.2, 2.8), (4.9, 3.1), (5.5, 2.6)]
        >>> materials = ['Material A', 'Material B', 'Material C']
        >>> energy_band_alignment_diagram(energies, materials, limit=8.0, width=0.8,
                                    cols=['#74356C', '#efce19'], textsize=18,
                                    arrowhead=0.5, outfile='BandAlignment',
                                    references=[(3.0, 'Reference 1'), (4.0, 'Reference 2')],
                                    edge='black')
    """
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
        ax1.bar(i,-limit, width, edgecolor=None)
        ax1.bar(i,-energies[i][1], width, color='w', edgecolor=None)

    ## Reset the colours back to the start and plot the EA
    ax1.set_color_cycle(cols)
    for i in ind:
        ax1.bar(i,-energies[i][0], width, edgecolor=None,alpha=0.8)

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
        ax1.arrow(i-0.25,-en[0],0, en[0]-arrowhead, width=0.005,
                  head_length=arrowhead,head_width=0.07, fc='black',ec='None')
        ax1.arrow(i-0.25,0,0, -en[1]+arrowhead, width=0.005,
                  head_length=arrowhead,head_width=0.07, fc='black',ec='None')
        ax1.arrow(i-0.25,0,0, -en[0]+arrowhead, width=0.005,
                  head_length=arrowhead,head_width=0.07, fc='black',ec='None')
        loc_ip = -(en[0] + en[1]) / 2
        ax1.text(i-os1,loc_ip,"IP  %3.1f"%en[1],fontsize=textsize)

        loc_ea = -en[0] / 2
        ax1.text(i-os2,loc_ea,"EA %3.1f"%en[0],fontsize=textsize)

        ax1.minorticks_on()
        # Don't show minor ticks on x-axis
        ax1.tick_params(axis='x',which='minor',bottom='off')
        ax2.minorticks_on()

    for ref in references:
        ax1.hlines(-ref[1], -0.5, len(energies) - 0.5,
                   linestyles='--', colors='r')
        ax1.text(len(energies) - 0.45, -ref[1] - 0.1, ref[0],
                 fontsize=textsize, color='r')

    fig.savefig('%s.eps'%outfile,bbox_inches='tight')
    fig.savefig('%s.png'%outfile,bbox_inches='tight')
    plt.show()
    print("Figure saved as %s.eps and %s.png"%(outfile, outfile))
    plt.close(fig)
