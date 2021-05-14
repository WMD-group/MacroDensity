#! /usr/bin/env python
#------------------------------------------------------------------------------
'''
def vasp_params():
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    return
'''
#------------------------------------------------------------------------------

def bulk_interstitial_alignment(interstices,outcar="OUTCAR",locpot="LOCPOT",cube_size=[1,1,1]):
    '''
    Alignment of the band edges with the interstitial "vacuum" potential.

    Inputs:
    intersices = Positions of the pores/interstices ([[interstice1],[interstice2],...])
    outcar = a VASP OUTCAR (str)
    locpot = a VASP LOCPOT (str)
    cube_size = a cube defined by LOCPOT FFT mesh points ([int,int,int] - optional)

    Output: Aligned Valence Band, Aligned Conduction Band
    '''
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, volume_average
    from macrodensity.vasp_tools import get_band_extrema

    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(locpot,quiet=True)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)

    band_extrema = get_band_extrema(outcar)
    VB_eigenvalue = band_extrema[0]
    CB_eigenvalue = band_extrema[1]
    ## Calculation of reference state from LOCPOT
    ## Avoiding list comprehensions to minimise the computational overhead.
    interstitial_potentials = []
    interstitial_variances = []
    for interstice in interstices:
        #origin, cube, grid, nx, ny, nz, travelled=[0, 0, 0]
        locpot_extract = volume_average(origin=interstice,cube=cube_size,grid=grid_pot,nx=NGX,ny=NGY,nz=NGZ)
        interstitial_potentials.append(locpot_extract[0])
        interstitial_variances.append(locpot_extract[1])
    ## Calculating the referenced band energies, then rounding the output
    ## REMINDER: Add variance output for further (verbose) analysis options
    sum_interstitial_potential = 0
    for ele in interstitial_potentials:
        sum_interstitial_potential += ele
    average_interstitial_potential = sum_interstitial_potential/len(interstitial_potentials)
    VB_aligned = round(VB_eigenvalue - average_interstitial_potential,2)
    CB_aligned = round(CB_eigenvalue - average_interstitial_potential,2)
    ## Data formatting and output
    print("Reading band edges from file: "+str(outcar))
    print("Reading potential from file: "+str(locpot))
    print("Interstital variances: "+str(interstitial_variances))
    print("VB_aligned      CB_aligned")
    print("--------------------------------")
    print(VB_aligned,"         ",CB_aligned)
    return VB_aligned, CB_aligned
#------------------------------------------------------------------------------

def plot_active_space(a_point,b_point,c_point,tolerance=1E-4,input_file='LOCPOT'):
    '''
    Input section (define the plane with 3 points)
    a_point = [0, 0, 0]
    b_point = [1, 0, 1]
    c_point = [0, 1, 0]
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    import csv
    ##from itertools import izip
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid, planar_average, macroscopic_average,cube_potential
    from macrodensity.beta_tools import create_plotting_mesh,points_2_plane
    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    cutoff_varience = tolerance
    ## Get the gradiens (Field), if required.
    ## Comment out if not required, due to compuational expense.
    grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    #------------------------------------------------------------------
    ## Get the equation for the plane
    ## This is the section for plotting on a user defined plane;
    ## uncomment commands if this is the option that you want.
    ##------------------------------------------------------------------
    ## Convert the fractional points to grid points on the density surface
    a = numbers_2_grid(a_point,NGX,NGY,NGZ)
    b = numbers_2_grid(b_point,NGX,NGY,NGZ)
    c = numbers_2_grid(c_point,NGX,NGY,NGZ)
    plane_coeff = points_2_plane(a,b,c)
    ## Get the gradients
    XY = np.multiply(grad_x,grad_y)
    grad_mag = np.multiply(XY,grad_z)
    ## Create the plane
    xx,yy,grd = create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_x)
    ## Plot the surface
    plt.contourf(xx,yy,grd,V)
    plt.show()
    ##------------------------------------------------------------------
    ## Plotting a planar average (Field/potential) throughout the material
    ##------------------------------------------------------------------
    ## FIELDS
    planar = planar_average(grad_x,NGX,NGY,NGZ)
    ## POTENTIAL
    planar = planar_average(grid_pot,NGX,NGY,NGZ)
    ## MACROSCOPIC AVERAGE
    macro  = macroscopic_average(planar,4.80,resolution_z)
    plt.plot(planar)
    plt.plot(macro)
    plt.savefig('Planar.eps')
    plt.show()
    ##------------------------------------------------------------------
    # Getting the average potential in a single cube of arbitrary size
    ##------------------------------------------------------------------
    ## cube defines the size of the cube in units of mesh points (NGX/Y/Z)
    cube = [2,2,2]
    ## origin defines the bottom left point of the cube the "0,0,0" point in fractional coordinates
    origin = [0,0,0]
    ## travelled; do not alter this variable
    travelled = [0,0,0]
    ## Uncomment the lines below to do the business
    vacuum = []
    non_vacuum = []
    for i in range(0,NGX,cube[0]):
        print(float(i)/NGX)
        for j in range(0,NGY,cube[1]):
            for k in range(0,NGZ,cube[2]):
                origin = [float(i)/NGX,float(j)/NGY,float(k)/NGZ]
                cube_potential, cube_var = cube_potential(origin,travelled,cube,grid_pot,NGX,NGY,NGZ)
            if cube_var <= cutoff_varience:
                vacuum.append(origin)
            else:
                non_vacuum.append(origin)
    print("Number of vacuum cubes: ", len(vacuum))
    print("Number of non-vacuum cubes: ", len(non_vacuum))
    print("Percentage of vacuum cubes: ",(float(len(vacuum))/(float(len(vacuum))+float(len(non_vacuum)))*100.))
    print("Percentage of non-vacuum cubes: ",(float(len(non_vacuum))/(float(len(vacuum))+float(len(non_vacuum)))*100.))
#------------------------------------------------------------------------------

def plot_field_at_point(a_point,b_point,c_point,input_file='LOCPOT'):
    '''
    Input section (define the plane with 3 points, fractional coordinates)
    a_point = [0, 0, 0]
    b_point = [1, 0, 1]
    c_point = [0, 1, 0]
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import colors,cm #colour maps; so I can specify cube helix
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid, planar_average, macroscopic_average,cube_potential
    from macrodensity.beta_tools import create_plotting_mesh,points_2_plane
    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    ## Get the gradiens (Field), if required.
    ## Comment out if not required, due to compuational expense.
    print("Calculating gradients (Electic field, E=-Grad.V )...")
    grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    #------------------------------------------------------------------
    ## Get the equation for the plane
    ## This is the section for plotting on a user defined plane;
    ## uncomment commands if this is the option that you want.
    ##------------------------------------------------------------------
    ## Convert the fractional points to grid points on the density surface
    a = numbers_2_grid(a_point,NGX,NGY,NGZ)
    b = numbers_2_grid(b_point,NGX,NGY,NGZ)
    c = numbers_2_grid(c_point,NGX,NGY,NGZ)
    plane_coeff = points_2_plane(a,b,c)
    ## Calculate magnitude of gradient.
    # Should be able to use numpy.linalg.norm for this, but the Python array indices are causing me grief
    X2 = np.multiply(grad_x,grad_x)
    Y2 = np.multiply(grad_y,grad_y)
    Z2 = np.multiply(grad_z,grad_z)
    grad_mag = np.sqrt(np.add(X2,Y2,Z2))
    # This was my, non working, attempt to use the built in function.
    #grad_mag=np.linalg.norm( [grad_y,grad_y,grad_z], axis=3)

    ## This function in Macrodensity averages Efield ACROSS Z for Slab calculations
    #xx,yy,grd =  pot.create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_mag) #AVG over full volume

    # Here we construct the same xx,yy,grd variables with a SLICE, forming a plane in XY at particular ZSLICE
    xx, yy = np.mgrid[0:NGX,0:NGY]
    ZSLICE= NGZ/2 # Chosses where in the Z axis the XY slice is cut through
    # Slice of magnitude of electric field, for contour plotting
    grd=grad_mag[:,:,ZSLICE]
    # Slices of x and y components for arrow plotting
    grad_x_slice=grad_x[:,:,ZSLICE]
    grad_y_slice=grad_y[:,:,ZSLICE]
    # OK, that's all our data

    # This code tiles the data to (2,2) to re-expand unit cell to a 2x2 supercell in XY
    xx,yy=np.mgrid[0:2*NGX,0:2*NGY]
    grd=np.tile(grd, (2,2))
    grad_x_slice=np.tile(grad_x_slice, (2,2))
    grad_y_slice=np.tile(grad_y_slice, (2,2))
    # End of tiling code

    ## Contours of the above sliced data
    plt.contour(xx,yy,grd,6,cmap=cm.cubehelix)

    # Also generate a set of Efield arrows ('quiver') for this data.
    # Specifying the drawing parameters is quite frustrating - they are very brittle + poorly documented.
    plt.quiver(xx,yy, grad_x_slice, grad_y_slice,
            color='grey',
            units='dots', width=1, headwidth=3, headlength=4
            ) #,
    #        units='xy', scale=10., zorder=3, color='blue',
    #        width=0.007, headwidth=3., headlength=4.)

    plt.axis('equal') #force square aspect ratio; this assuming X and Y are equal.
    plt.show()
#------------------------------------------------------------------------------

def plot_gulp_potential(lattice_vector,input_file='gulp.out',output_file='planar.dat',new_resolution = 3000):
    '''
    lattice_vector = 3.0
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d
    from macrodensity.density_tools import matrix_2_abc, planar_average, macroscopic_average,density_2_grid_gulp,read_gulp_potential
    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    pot, NGX, NGY, NGZ, Lattice = read_gulp_potential(input_file)
    vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot = density_2_grid_gulp(pot, NGX, NGY, NGZ)
    #------------------------------------------------------------------
    ## POTENTIAL PLANAR AVERAGE
    planar = planar_average(grid_pot, NGX, NGY, NGZ)
    np.savetxt(output_file, planar)
    ## MACROSCOPIC AVERAGE
    new_abscissa = np.linspace(0, NGZ - 1, new_resolution)
    f = interp1d(range(NGZ), planar, kind='cubic')
    interpolated_potential = [f(i) for i in new_abscissa]
    macro  = macroscopic_average(planar, lattice_vector, vector_c/new_resolution)
    ## PLOT
    plt.plot(interpolated_potential)
    plt.show()
    plt.plot(planar)
    plt.show()
#------------------------------------------------------------------------------

def plot_on_site_potential(species,sample_cube,potential_file='LOCPOT',coordinate_file='POSCAR'):
    '''
    potential_file = 'LOCPOT' # The file with VASP output for potential
    coordinate_file = 'POSCAR' # The coordinates file NOTE NOTE This must be in vasp 4 format
    species = "O"  # The species whose on-site potential you are interested in
    sample_cube = [5,5,5] # The size of the sampling cube in units of mesh points (NGX/Y/Z)
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    import csv
    #from itertools import izip
    import ase                # Only add this if want to read in coordinates
    from ase.io import write  # Only add this if want to read in coordinates
    from ase.io import vasp   # Only add this if want to read in coordinates
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid, planar_average, macroscopic_average,volume_average
    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(potential_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    ## Get the gradiens (Field), if required.
    ## Comment out if not required, due to compuational expense.
    grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    ##------------------------------------------------------------------
    ## Getting the potentials for a group of atoms, in this case the Os
    ##------------------------------------------------------------------
    coords = ase.io.vasp.read_vasp(coordinate_file)
    scaled_coords = coords.get_scaled_positions()
    symbols = coords.get_chemical_symbols()
    ox_coords = []

    for i, atom in enumerate(coords):
        if symbols[i] == species:
            ox_coords.append(scaled_coords[i])
    grid_position = np.zeros(shape=(3))
    potentials_list = []
    i = 0
    num_bins = 20
    for coord in ox_coords:
        i = i + 1
        grid_position[0] = coord[0]
        grid_position[1] = coord[1]
        grid_position[2] = coord[2]
        cube = sample_cube    # The size of the cube x,y,z in units of grid resolution.
        origin = [grid_position[0]-2,grid_position[1]-2,grid_position[2]-1]
        travelled = [0,0,0] # Should be left as it is.
        cube_potential, cube_var = volume_average(origin,cube,grid_pot,NGX,NGY,NGZ)
        potentials_list.append(cube_potential)
    n, bins, patches = plt.hist(potentials_list, num_bins, facecolor='#6400E1', alpha=0.5) ##normed=100
    plt.xlabel('Hartree potential (V)',fontsize = 22)
    plt.ylabel('% of centres',fontsize = 22)
    plt.savefig('Potentials.png',dpi=300)
    plt.show()
#------------------------------------------------------------------------------

def plot_planar_average(lattice_vector,input_file="LOCPOT",output_file="Planar.out"):
    import numpy as np
    import matplotlib.pyplot as plt
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, planar_average, macroscopic_average
    # Get the potential
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    ## POTENTIAL
    planar = planar_average(grid_pot,NGX,NGY,NGZ)
    ## MACROSCOPIC AVERAGE
    macro = macroscopic_average(planar,lattice_vector,resolution_z)
    plt.plot(planar)
    plt.plot(macro)
    plt.savefig('Planar.eps')
    plt.show()
    np.savetxt(output_file,planar)
#------------------------------------------------------------------------------

def plot_planar_cube(input_file,lattice_vector,output_file = 'planar.dat'):
    '''
    input_file = 'cube_002_total_density.cube'
    lattice_vector = 4.75
    output_file = 'planar.dat'
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    import ase.io.cube
    from macrodensity.density_tools import planar_average, macroscopic_average
    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    potential, atoms = ase.io.cube.read_cube(input_file,read_data=True)
    vector_a = np.linalg.norm(atoms.cell[1])
    vector_b = np.linalg.norm(atoms.cell[1])
    vector_c = np.linalg.norm(atoms.cell[2])
    NGX = len(potential)
    NGY = len(potential[0])
    NGZ = len(potential[0][0])
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    print(NGX,NGY,NGZ)
    #------------------------------------------------------------------
    ## POTENTIAL
    planar = planar_average(potential,NGX,NGY,NGZ)
    ## MACROSCOPIC AVERAGE
    macro  = macroscopic_average(planar,lattice_vector,resolution_z)
    plt.plot(planar)
    plt.plot(macro)
    plt.savefig('Planar.eps')
    #plt.show()
    np.savetxt(output_file,planar)
#------------------------------------------------------------------------------

def plot_plane_field(a_point,b_point,c_point,input_file='LOCPOT'):
    '''
    Input section (define the plane with 3 points, fractional coordinates)
    a_point = [0, 0, 0]
    b_point = [1, 0, 1]
    c_point = [0, 1, 0]
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid
    from macrodensity.beta_tools import create_plotting_mesh,points_2_plane
    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    ## Get the gradiens (Field), if required.
    ## Comment out if not required, due to compuational expense.
    grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    #------------------------------------------------------------------
    ## Get the equation for the plane
    ## This is the section for plotting on a user defined plane;
    ## uncomment commands if this is the option that you want.
    ##------------------------------------------------------------------
    ## Convert the fractional points to grid points on the density surface
    a = numbers_2_grid(a_point,NGX,NGY,NGZ)
    b = numbers_2_grid(b_point,NGX,NGY,NGZ)
    c = numbers_2_grid(c_point,NGX,NGY,NGZ)
    plane_coeff = points_2_plane(a,b,c)
    ## Get the gradients
    XY = np.multiply(grad_x,grad_y)
    grad_mag = np.multiply(XY,grad_z)
    ## Create the plane
    xx,yy,grd = create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_x)
    ## Plot the surface
    plt.contour(xx,yy,grd,1)
    plt.show()
#------------------------------------------------------------------------------
def moving_cube(cube=[1,1,1],vector=[1,1,1],origin=[0,0,0],magnitude = 280,input_file='LOCPOT'):
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, vector_2_abscissa, travelling_volume_average
    import matplotlib.pyplot as plt
    '''
    ## cube defines the size of the cube in units of mesh points (NGX/Y/Z)
    ## vector is the vector you wish to travel along
    ## origin defines the origin of the line in units of mesh points (NGX/Y/Z)
    ## magnitude defines the length of the line, in units of mesh points (NGX/Y/Z)
    '''
    #------------------------------------------------------------------
    # Get the potential
    # This section should not be altered
    #------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)

    ##------------------------------------------------------------------
    ## Plotting the average in a moving cube along a vector
    ##------------------------------------------------------------------

    ## IF YOU WANT TO PLOT THE POTENTIAL:
    cubes_potential = travelling_volume_average(grid_pot,cube,origin,vector,NGX,NGY,NGZ,magnitude)
    abscissa = vector_2_abscissa(vector,magnitude,resolution_x,resolution_y,resolution_z)
    plt.plot(abscissa, cubes_potential)
    plt.xlabel("$z (\AA)$")
    plt.ylabel("Potential (eV)")
    plt.savefig('moving_cube.png')
#------------------------------------------------------------------------------
def spherical_average(cube_size,cube_origin,input_file='LOCPOT'):
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, volume_average
    """ Calculates the Spherical Average
    Inputs:
    input_file = i.e. 'LOCPOT'
    cube_size = e.g. [2,2,2] This size is in units of mesh points (NGX/Y/Z)
    cube_origin = e.g. [0,0,0] Defines the bottom left point of the cube the "0,0,0" point in fractional coordinates
    Output: cube_potential, cube_variance
    """
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)

    cube = cube_size
    origin = cube_origin
    ## travelled; do not alter this variable
    travelled = [0,0,0]
    cube_pot, cube_var = volume_average(origin=cube_origin,cube=cube_size,grid=grid_pot,nx=NGX,ny=NGY,nz=NGZ,travelled=[0,0,0])
    print("Potential            Variance")
    print("--------------------------------")
    print(cube_pot,"   ", cube_var)
    return cube_pot, cube_var
#------------------------------------------------------------------------------
