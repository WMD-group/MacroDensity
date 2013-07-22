#!/bin/bash
import numpy
import matplotlib.pyplot as plt
from pylab import *
import math
from matplotlib.patches import Polygon
def distance(a, b):
    # distance between points
    return numpy.sqrt(numpy.sum((a - b)**2))

def calc_sphere_field(NGX,NGY,NGZ,lattice,Potential,centre_sphere,radius,cutoff_field,f):
    """Calculate the electric field across the sphere in a number of directions"""
    x_resolution = lattice[0,0]/NGX
    y_resolution = lattice[0,0]/NGY
    z_resolution = lattice[0,0]/NGZ
    centre =  np.zeros(shape=(3))
    rho_a =  np.zeros(shape=(3))
    rho_a_prime =  np.zeros(shape=(3))
    centre[0] = int(round(centroid[0]*NGX))
    centre[1] = int(round(centroid[1]*NGY))
    centre[2] = int(round(centroid[2]*NGZ))

    for x in range (0,int(round(radius/x_resolution))):
        for y in range (0, int(round((radius - (x*x_resolution)**3)/y_resolution))):
            z = int(round((radius - x*(x_resolution)**3 - (y*y_resolution)**3)/z_resolution))
	    rho_a = Potential[centre[0]+x-int((centre[0]+x)/NGX)*NGX,centre[1]+y-int((centre[1]+y)/NGY)*NGY,centre[2]+z-int((centre[2]+z)/NGZ)*NGZ]
	    rho_a_prime = Potential[centre[0]-x+int((centre[0]-x)/NGX)*NGX,centre[1]-y+int((centre[1]-y)/NGY)*NGY,centre[2]-z+int((centre[2]-z)/NGZ)*NGZ]
	    Field = linalg.norm(rho_a - rho_a_prime)/(2*radius)

	    if Field > cutoff_field:
 	        f.write("Warning, field across sphere exceeds cutoff :")
	        f.write(str(Field))
		f.write('\n')
	        print("Warning, field across sphere exceeds cutoff :", Field)

def macro_av(NGX,NGY,NGZ,Plane_Potential_New):
    # Macroscopic Averaging
    #
    #
    Macro_Potential = numpy.zeros(shape=(NGZ))
    Macro_Potential_B = numpy.zeros(shape=(NGZ))
    Period = float(raw_input('What periodicity do you want for macroscopic averaging?  '))
    Period_Points = int(Period/(lattice[2,2]/NGZ))
    for i in range (0, NGZ):
     k = 0
     Macro_Potential[i] = 0
     for j in range(i-Period_Points/2,i+Period_Points/2+1):
      if j <  NGZ and j >= 0:
       Macro_Potential[i] = Macro_Potential[i] + Plane_Potential_New[j]
      if j >=  NGZ:
       Macro_Potential[i] = Macro_Potential[i] + Plane_Potential_New[j-NGZ]
      if j <  0:
       Macro_Potential[i] = Macro_Potential[i] + Plane_Potential_New[j+NGZ]
     Macro_Potential[i] = Macro_Potential[i]/Period_Points
    for i in range (0, NGZ):
     k = 0
     Macro_Potential_B[i] = 0
     for j in range(i-Period_Points/2,i+Period_Points/2+1):
      if j <  NGZ and j >= 0:
       Macro_Potential_B[i] = Macro_Potential_B[i] + Macro_Potential[j]
      if j >=  NGZ:
       Macro_Potential_B[i] = Macro_Potential_B[i] + Macro_Potential[j-NGZ]
      if j <  0:
       Macro_Potential_B[i] = Macro_Potential_B[i] + Macro_Potential[j+NGZ]
     Macro_Potential_B[i] = Macro_Potential_B[i]/Period_Points
    
    return Macro_Potential_B

# Sphere_Potential = spher_av(NGX,NGY,NGZ,Potential_grid,axis,centroid,lattice)
def spher_av(NGX,NGY,NGZ,Potential,axis,centroid,lattice):
    # Spherical averaging
    radius = float(raw_input("What radius would you like for spherical averaging? "))
    centre_sphere = numpy.zeros(shape=(3))
    point = numpy.zeros(shape=(3))
    dpoint = numpy.zeros(shape=(3))
    sphere_pot = numpy.zeros(shape=(NGZ))
    spherical_average=numpy.zeros(shape=(NGZ,2))
    for z in range (0, NGZ):
     sphere_pot_list = []
     centre_sphere[0] = lattice[0,0] * centroid[0]
     centre_sphere[1] = lattice[1,1] * centroid[1]
     centre_sphere[2] = lattice[2,2]/NGZ * z
# If the centre is more than radius away from the edge, restrict the search
     if z*lattice[2,2]/NGZ > radius and z*lattice[2,2]/NGZ < lattice[2,2] - radius:
      for i in range (0, NGX):
       for j in range (0, NGY):
        for k in range (z-int(radius*NGZ/lattice[2,2]), z+int(radius*NGZ/lattice[2,2])):
         point[0] = lattice[0,0]/NGX * i
         point[1] = lattice[1,1]/NGY * j
         point[2] = lattice[2,2]/NGZ * k 
         separation = distance(point,centre_sphere) 
         if separation <= radius:
          sphere_pot_list.append(Potential[i,j,k])
     else:
      for i in range (0, NGX):
       for j in range (0, NGY):
        for k in range (0, NGZ):
	 point[0] = lattice[0,0]/NGX * i
	 point[1] = lattice[1,1]/NGY * j
	 point[2] = lattice[2,2]/NGZ * k
# Minimum image convention
         dpoint[0] = point[0] - centre_sphere[0]
         dpoint[0] = dpoint[0] - round(dpoint[0]/lattice[0,0])*lattice[0,0]
         dpoint[1] = point[1] - centre_sphere[1]
         dpoint[1] = dpoint[1] - round(dpoint[1]/lattice[1,1])*lattice[1,1]
         dpoint[2] = point[2] - centre_sphere[2]
         dpoint[2] = dpoint[2] - round(dpoint[2]/lattice[2,2])*lattice[2,2]
 	 separation = numpy.sqrt(dpoint[0]**2 + dpoint[1]**2 + dpoint[2]**2)
         if separation <= radius:
          sphere_pot_list.append(Potential[i,j,k])
    
     print(z, numpy.mean(sphere_pot_list), numpy.var(sphere_pot_list))
     spherical_average[z,0] = numpy.mean(sphere_pot_list)
     spherical_average[z,1] = numpy.var(sphere_pot_list)
    return spherical_average

def point_sphere(NGX,NGY,NGZ,Potential,axis,centroid,lattice,f):
    """Calculates the spherical average, only at given points""" 
    radius = float(raw_input("What radius would you like for spherical averaging? "))
    cutoff_field = float(raw_input("What is the cutoff value for electric field across the sphere? (V/A) "))
    centre_sphere = numpy.zeros(shape=(3))
    point = numpy.zeros(shape=(3))
    dpoint = numpy.zeros(shape=(3))
    sphere_pot = numpy.zeros(shape=(NGZ))
    spherical_average=numpy.zeros(shape=(2))
    sphere_pot_list = []
    centre_sphere[0] = lattice[0,0] * centroid[0]
    centre_sphere[1] = lattice[1,1] * centroid[1]
    centre_sphere[2] = lattice[2,2] * centroid[2]
# If the centre is more than radius away from the edge, restrict the search
    if centre_sphere[2] > radius and centre_sphere[2] < lattice[2,2] - radius:
     for i in range (0, NGX):
      for j in range (0, NGY):
       for k in range (int((centre_sphere[2]-radius)*NGZ/lattice[2,2]), int((centre_sphere[2]+radius)*NGZ/lattice[2,2])):
        point[0] = lattice[0,0]/NGX * i
        point[1] = lattice[1,1]/NGY * j
        point[2] = lattice[2,2]/NGZ * k
        separation = distance(point,centre_sphere)
        if separation <= radius:
         sphere_pot_list.append(Potential[i,j,k])

    else:
     for i in range (0, NGX):
      for j in range (0, NGY):
       for k in range (0, NGZ):
         point[0] = lattice[0,0]/NGX * i
         point[1] = lattice[1,1]/NGY * j
         point[2] = lattice[2,2]/NGZ * k
# Minimum image convention
         dpoint[0] = point[0] - centre_sphere[0]
         dpoint[0] = dpoint[0] - round(dpoint[0]/lattice[0,0])*lattice[0,0]
         dpoint[1] = point[1] - centre_sphere[1]
         dpoint[1] = dpoint[1] - round(dpoint[1]/lattice[1,1])*lattice[1,1]
         dpoint[2] = point[2] - centre_sphere[2]
         dpoint[2] = dpoint[2] - round(dpoint[2]/lattice[2,2])*lattice[2,2]
         separation = numpy.sqrt(dpoint[0]**2 + dpoint[1]**2 + dpoint[2]**2)
         if separation <= radius:
          sphere_pot_list.append(Potential[i,j,k])

    field = calc_sphere_field(NGX,NGY,NGZ,lattice,Potential,centre_sphere,radius,cutoff_field,f)
    spherical_average[0] = numpy.mean(sphere_pot_list)
    spherical_average[1] = numpy.var(sphere_pot_list)

    return spherical_average


def list_2_matrix(Potential,NGX,NGY,NGZ):
    # Convert the linear list of numbers to "cartesinan" grid
    i = 0
    j = 0
    k = 0
    Potential_grid = numpy.zeros(shape=(NGX,NGY,NGZ))
    for i in range (0, NGZ):
     for j in range (0, NGY):
      for k in range (0, NGX):
       Potential_grid[k,j,i] = Potential[k+(j)*NGX+(i)*NGY*NGX]

    return Potential_grid

def planar_av(NGX,NGY,NGZ,Potential):
    # Planar average potential
    i = 0
    j = 0
    k = 0
    Plane_Potential = numpy.zeros(shape=(NGZ))
    for i in range (0, NGZ):
     for j in range (0, NGY):
      for k in range (0, NGX):
       Plane_Potential[i]=Plane_Potential[i]+Potential[k+(j-1)*NGX+(i-1)*NGY*NGX]
       k = k + 1
      j = j + 1
     Plane_Potential[i] = Plane_Potential[i]/(NGX*NGY)
     i = i + 1
    return Plane_Potential

def shift_plot(Plane_Potential,NGZ,bin_shift):
    # Now move the potentials by the required shift
    #
    #
    Plane_Potential_New = numpy.zeros(shape=(NGZ))
    i = 0
    
    for i in range (0, NGZ):
     if i + bin_shift >= 0 and i + bin_shift < NGZ:
      Plane_Potential_New[i+bin_shift] = Plane_Potential[i]
     if i + bin_shift < 0:
      Plane_Potential_New[i+bin_shift+NGZ] = Plane_Potential[i]
     if i + bin_shift >= NGZ:
      Plane_Potential_New[i+bin_shift-NGZ] = Plane_Potential[i]
    

    return Plane_Potential_New

# Open the input
f = open('LOCPOT',"r")
lines = f.readlines()
f.close()
# Split the LOCPOT into relevant sections

lattice = numpy.zeros(shape=(3,3))
NGX = 0
NGY = 0
NGZ = 0
i=0
# Get Header information
for line in lines:
 inp = line.split()
 if inp == []:
  continue
 if len(inp) > 0:
  i = i+1
 if i >= 3 and i < 6:
  lattice[i-3,:]=inp[:]
 if i == 6:
  num_species=len(inp)
  species=inp
 if i == 7:
  num_type=inp
  j = 0
  while (j < num_species):
   num_type[j-1] = int(num_type[j-1])
   j = j + 1
  num_atoms=sum(num_type)
 if i == 8:
  coord_type = inp
 
Surface = sum(numpy.cross(lattice[0,:],lattice[1,:]))
# Restart reading to get the coordinates...it's just easier this way!
i=0
Coordinates = numpy.zeros(shape=(num_atoms,3))
for line in lines:
 inp = line.split()
 if len(inp) > 0:
  i = i + 1
 if i >= 9 and i <= num_atoms+8 and len(inp) > 0:
  Coordinates[i-9,0] = float(inp[0])
  Coordinates[i-9,1] = float(inp[1])
  Coordinates[i-9,2] = float(inp[2])
# Now get the info about the charge grid
i = 0
for line in lines:
 inp = line.split()
 if len(inp) > 0:
  i = i + 1
 if i == num_atoms + 9:
  NGX = int(inp[0])
  NGY = int(inp[1])
  NGZ = int(inp[2])
  k = 0
  Potential = numpy.zeros(shape=(NGX*NGY*NGZ))
# Read in the potential data
 if i > num_atoms + 9 and i < num_atoms + 9 + NGX*NGY*NGZ/5:
  Potential[k]   = inp[0]
  Potential[k+1] = inp[1]
  Potential[k+2] = inp[2]
  Potential[k+3] = inp[3]
  Potential[k+4] = inp[4]
  k = k + 5
  if math.fmod(k,100000) == 0:
   print "Reading potetial, at point", k

# Start processing the data
average_type=raw_input("Which kind of average would you like? (P)lanar/(S)pherical/(Po)int spherical average ")
if average_type != "Po":
 axis = raw_input("Which axis do you want to plot along? (X/Y/Z)")

# Default axis is set to Z, if X or Y are required we need to re-jig things.
#Since the code was all written for Z I just rename the different parameters.
#Its a dirty little fix, and I like it!
 if axis == 'Y':
  tmp = NGZ ; NGZ = NGY ; NGY = NGX ; NGX = tmp
  tmp = lattice[2,2] ; lattice[2,2] = lattice[1,1] ; lattice[1,1] = lattice[0,0]
  lattice[0,0] = tmp
 elif axis == 'X':
  tmp = NGZ ; NGZ = NGX ; NGX = NGY ; NGY = tmp
  tmp = lattice[2,2] ; lattice[2,2] = lattice[0,0] ; lattice[0,0] = lattice[1,1]
  lattice[1,1] = tmp

#3 Convert the data to a more friendly format

if average_type != "P":
 Potential_grid = list_2_matrix(Potential,NGX,NGY,NGZ)

if average_type == 'P':
 f = open('Planar_Av.dat','w')
 f.write('Planar Average of Potential \n')
 Plane_Potential = planar_av(NGX,NGY,NGZ,Potential)
 # Section for re-centering the plot, so they can be consistent
 Centre = float(raw_input('Where do you want the plot centred?  '))
 # How many bins do we need to shift by?
 bin_shift =int((lattice[2,2]/2 - Centre) / (lattice[2,2] / NGZ))
 # Now move the potentials by the required shift
 Plane_Potential_New = shift_plot(Plane_Potential,NGZ,bin_shift)
 np.savetxt(f,Plane_Potential_New)
 # Now get the macroscopic average
 f.write('Macroscopic Average of Potential \n')
 Macro_Potential = macro_av(NGX,NGY,NGZ,Plane_Potential_New)
 np.savetxt(f,Macro_Potential)
 Vacuum_potential = Macro_Potential[0]
 Centre_potential = Macro_Potential[NGZ/2]
 f.write("Average bulk potential at centre of slab: ")
 f.write(Centre_potential)
 f.write('\n')
 f.write("Average potential of vacuum : vacuum_potential: ")
 f.write(Vacuum_potential)
 f.write('\n')
elif average_type == 'S':
 f = open('Spherical_Av.dat','w')
 f.write('Spherical Average of Potential \n')
 spherical_average=numpy.zeros(shape=(NGZ,2))
 spherical_av_potential=numpy.zeros(shape=(NGZ))
 centroid = numpy.zeros(shape=(2))
 centroid[0] = raw_input("a value on the plane to centre the sphere (0 - 1) ")
 centroid[1] = raw_input("b value on the plane to centre the sphere (0 - 1) ")
 spherical_average = spher_av(NGX,NGY,NGZ,Potential_grid,axis,centroid,lattice)
 np.savetxt(f,spher_av)
 for i in range (0,NGZ):
  spherical_av_potential[i] = spherical_average[i,0]
elif average_type == 'Po':
 f = open('Spherical_Av.dat','w')
 f.write('Spherical Average of Potential \n')
 number_points = raw_input("How many points would you like to consider? ")
 for i in range (0,int(number_points)):
  spherical_average=numpy.zeros(shape=(NGZ,2))
  spherical_av_potential=numpy.zeros(shape=(NGZ))
  centroid = numpy.zeros(shape=(3))
  centroid[0] = raw_input("a lattice value to centre the sphere (0 - 1) ")
  centroid[1] = raw_input("b lattice value to centre the sphere (0 - 1) ")
  centroid[2] = raw_input("c lattice value to centre the sphere (0 - 1) ")
  f.write('Spherical Average of Potential at point \n')
  np.savetxt(f,centroid)
  f.write('\n')
  spherical_average = point_sphere(NGX,NGY,NGZ,Potential_grid,axis,centroid,lattice,f)
  f.write("   Average    Variance \n")
  np.savetxt(f,spherical_average)
  print("   Average    Variance")
  print(spherical_average)
  i = i + 1
# Z axis values
ZAxis = numpy.zeros(shape=(NGZ))
i = 0
while (i < NGZ):
 ZAxis[i] = i*lattice[2,2]/NGZ
 i = i + 1


# Plot the planar average
#
#
# MATLAB style
if average_type != "Po":
 xticklines = getp(gca(), 'xticklines')
 yticklines = getp(gca(), 'yticklines')
 xgridlines = getp(gca(), 'xgridlines')
 ygridlines = getp(gca(), 'ygridlines')
 xticklabels = getp(gca(), 'xticklabels')
 ygridlines = getp(gca(), 'ygridlines')
 xticklabels = getp(gca(), 'xticklabels')
 yticklabels = getp(gca(), 'yticklabels')
 
 setp(xticklines, 'linewidth', 3)
 setp(yticklines, 'linewidth', 3)
#setp(xgridlines, 'linestyle', '-')
#setp(ygridlines, 'linestyle', '-')
 setp(yticklabels, 'color', 'Black', fontsize='medium')
 setp(xticklabels, 'color', 'Black', fontsize='medium')


 plt.xlabel('$z \AA$ ',fontsize='large')
 plt.grid(True)
 if average_type == 'P':
  plt.ylabel('$V_{planar}(z)$  $eV \AA^{-2} $',fontsize='large')
  plt.plot(ZAxis,Macro_Potential,ZAxis, Plane_Potential_New)
 elif average_type == 'S':
  plt.ylabel('$V_{spherical}(z)$  $eV \AA^{-2} $',fontsize='large')
  plt.plot(ZAxis,spherical_av_potential)

 plt.show()
 print("Central bulk potential: ",Centre_potential)
 print("Vacuum potential: ",Vacuum_potential)
