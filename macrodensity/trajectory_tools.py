###############################################################################
# Copyright Keith Butler(2019)                                                #
#                                                                             #
# This file MacroDensity.trajectory_tools.py is free software: you can        #
# redistribute it and/or modify it under the terms of the GNU General Public  #
# License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version.                         #
# This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    #
# more details.                                                               #
# You should have received a copy of the GNU General Public License along with#
# this program. If not, see <http://www.gnu.org/licenses/>.                   #
#                                                                             #
###############################################################################

from __future__ import print_function
import numpy as np
import string

def str2list(rawstr):
    rawlist = rawstr.strip(string.whitespace).split(' ')
    # Remove space elements in list.
    cleanlist = [x for x in rawlist if x != ' ' and x != '']
    return cleanlist

def line2list(line, field=' ', dtype=float):
    "Convert text data in a line to data object list."
    strlist = line.strip().split(field)
    if type(dtype) != type:
        raise TypeError('Illegal dtype.')
    datalist = [dtype(i) for i in strlist if i != '']

    return datalist

class Trajectory:
    """
    The trajectory class creates an object from a molecular dynamics trajectory.
    """
    def __init__(self, source, file_format='vasp'):

        self.source = source 
        if file_format == 'vasp':
            self.read_vasp_trajectory()

    def read_vasp_trajectory(self):
        '''
        Reads in the XDATCAR format
        Args:
            filename : the file containing the trajectory
        Returns:
            nsteps
            symbols
            lattice
            trajectory
        '''

        with open(self.source, 'r') as f:
            data = f.readlines()
            filelength = len(data)
            data = []

        with open(self.source, 'r') as f:
            # Skip two lines
            _ = f.readline()
            _ = f.readline()

            # lattice basis
            self.lattice = []
            for i in range(3):
                basis = line2list(f.readline())
                self.lattice.append(basis)

            # atom info
            symbol_list = []
            self.symbols = []
            self.atom_types = str2list(f.readline())
            self.atoms_num = str2list(f.readline())
            for i, num in enumerate(self.atoms_num):
                symbol = [str(self.atom_types[i])] * int(num)
                symbol_list.append(symbol)
            
            self.symbols = [item for sublist in symbol_list for item in sublist]
            self.nsteps = int((filelength-6)/(len(self.symbols)+1))

            self.trajectory = []
            step = 0
            for step in range(self.nsteps):
                coords = []
                _ = f.readline()
                for i in range(len(self.symbols)):
                    line = f.readline()
                    if not line:
                        print('Trajectory not fully printed for step {}' % step)
                        break
                    else:
                        coords.append(line2list(line))
               
                self.trajectory.append(np.array(coords))

def poscars_from_trajectory(trajectory, number=1, random=True):
    '''
    Generate a set of input files from a trajectory
    Args:
        trajectory: a Trajectory object
        number: the number of poscars to create
        random: sample the trajectory at random? bool
    '''

    if random:
        import random
        indices = random.sample(range(0, trajectory.nsteps), number)
    else:
        indices = np.arange(0, trajectory.nsteps, int(trajectory.nsteps/number))

    for i in indices:
        with open('POSCAR_%s' % i, 'w+') as f:
            f.write('Created By Chance\n')
            f.write('1.0\n')
            for line in trajectory.lattice:
                f.write('  %10.8f   %10.8f   %10.8f \n' % (line[0], line[1], line[2]))
            for symbol in trajectory.atom_types:
                f.write('  %s  ' % symbol)
            f.write('\n')
            for num in trajectory.atoms_num:
                f.write('  %s  ' % num)
            f.write('\n')
            f.write(' Direct \n')
            for coord in trajectory.trajectory[i]:
                f.write('  %10.8f   %10.8f   %10.8f \n' % (coord[0], coord[1], coord[2]))






