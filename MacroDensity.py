#! /usr/bin/env python

import numpy
import matplotlib.pyplot as plt
from pylab import *
import math
from matplotlib.patches import Polygon

def distance(a, b):
    # distance between points
    return numpy.sqrt(numpy.sum((a - b)**2))

def lattice_vectors(latt):
    """Converts lattice matrix to a,b,c,alpha,beta,gamma"""
# TODO: add the angles part!
    out = np.zeros(shape=(3))
    for i in range(3):
     out[i] = sqrt(sum(latt[i,:]**2))
    return(out) 
   

def period(a,b,c,NGX,NGY,NGZ):
    """Wrap the points a,b,c so that they are withing the box NGXYZ"""
    abc = np.zeros(shape=(3))
    if a < 0:
	abc[0] = a + NGX
    elif a >= NGX:
	abc[0] = a - NGX
    else:
	abc[0] = a
    if b < 0:
	abc[1] = b + NGY
    elif b >= NGX:
	abc[1] = b - NGY
    else:
	abc[1] = b
    if c < 0:
	abc[2] = c + NGZ
    elif c >= NGZ:
	abc[2] = c - NGZ
    else:
	abc[2] = c
    return(abc)


def cent_potential(NGX,NGY,NGZ,P,centroid,latt):
    """Calculate the potential at the sphere centre"""
    centre = 0
    lat_a = np.sqrt(sum(latt[0,:]**2))
    lat_b = np.sqrt(sum(latt[1,:]**2))
    lat_c = np.sqrt(sum(latt[2,:]**2))
    x_res = lat_a/NGX
    y_res = lat_b/NGY
    z_res = lat_c/NGZ
    centre = float(P[centroid[0]/x_res*NGY*NGZ+centroid[1]/y_res*NGZ+centroid[2]/z_res,3]) 
    return centre


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
def spher_av(P,centroid,latt,NGZ):

    radius = float(raw_input("What radius would you like for spherical averaging? "))
    centre_sphere = numpy.zeros(shape=(3))
    point = numpy.zeros(shape=(3))
    dpoint = numpy.zeros(shape=(3))
    sphere_pot = numpy.zeros(shape=(NGZ))
    spherical_average=numpy.zeros(shape=(NGZ,2))
    lat_a = sum(latt[0,:]**2)
    lat_b = sum(latt[1,:]**2)
    lat_c = sum(latt[2,:]**2)
    for z in range (NGZ):
     sphere_pot_list = []
     centre_sphere[0] = centroid[0]
     centre_sphere[1] = centroid[1]
     centre_sphere[2] = latt[2,2]/NGZ * z
     for i in range (len(P)):
      point[0] = float(P[i,0]) - int((P[i,0]-centroid[0])/lat_a)*lat_a
      point[1] = float(P[i,1]) - int((P[i,1]-centroid[1])/lat_b)*lat_b
      point[2] = float(P[i,2]) - int((P[i,2]-centre_sphere[2])/lat_c)*lat_c
      separation = distance(point,centre_sphere)
      if separation <= radius:
#	  print(separation,radius,P[i,3])
          sphere_pot_list.append(P[i,3])


#     print(z, numpy.mean(sphere_pot_list), numpy.var(sphere_pot_list))
     spherical_average[z,0] = numpy.mean(sphere_pot_list)
     spherical_average[z,1] = numpy.var(sphere_pot_list)

    return spherical_average

def point_sphere(NGX,NGY,NGZ,P,axis,centroid,latt,f,r):
    """Calculates the spherical average, only at given points""" 
    lat = lattice_vectors(latt)
    print centroid
    Field=np.zeros(shape=(3))
#    cutoff_field = float(raw_input("What is the cutoff value for electric field across the sphere? (V/A) "))
# Set the limits of the box for searching, the +1 ensures that the sphere is totally enclosed
    xr = int(NGX*r/lat[0]) + 1; yr = int(NGY*r/lat[1]) + 1; zr = int(NGZ*r/lat[2]) + 1
    spherical_average=numpy.zeros(shape=(3))
    sphere_pot_list = []
    point = np.zeros(shape=(3))
    cent = np.zeros(shape=(3))
    c0 = int((centroid[0]*NGX/lat[0]))
    c1 = int((centroid[1]*NGY/lat[1]))
    c2 = int((centroid[2]*NGZ/lat[2]))
# If the centre is more than radius away from the edge, restrict the search
    for i in range (int(c0-xr),int(c0+xr)):
        for j in range (int(c1-yr),int(c1+yr)):
	    for k in range (int(c2-zr),int(c2+zr)):
                point[0] = i*lat[0]/NGX; point[1] = j*lat[1]/NGY; point[2] = k*lat[2]/NGZ
		cent[0] = c0*lat[0]/NGX; cent[1] = c1*lat[1]/NGY; cent[2] = c2*lat[2]/NGZ
		separation = distance(cent,point)
#		print(mesh[0]*NGY*NGZ+mesh[1]*NGZ+mesh[2])
#		print(point)
		if separation < r:
		    mesh = period(i,j,k,NGX,NGY,NGZ)
		    point[0]=P[mesh[0]*NGY*NGZ+mesh[1]*NGZ+mesh[2],0]
		    point[1]=P[mesh[0]*NGY*NGZ+mesh[1]*NGZ+mesh[2],1]
		    point[2]=P[mesh[0]*NGY*NGZ+mesh[1]*NGZ+mesh[2],2]
		    sphere_pot_list.append(P[mesh[0]*NGY*NGZ+mesh[1]*NGZ+mesh[2],3])

    Field = calc_field_tensors(c0,c1,c2,P,xr,yr,zr,NGX,NGY,NGZ)
    spherical_average[0] = P[c0*NGY*NGZ+c1*NGZ+c2,3]
    spherical_average[1] = numpy.mean(sphere_pot_list)
    spherical_average[2] = numpy.var(sphere_pot_list)
    print("E_xx: ",Field[0])
    print("E_yy: ",Field[1])
    print("E_zz: ",Field[2])

    return spherical_average


def  calc_field_tensors(c0,c1,c2,P,xr,yr,zr,NGX,NGY,NGZ):
    """Get the tensor of the Elctric field xx, yy, zz """
    F=np.zeros(shape=(3))
    mesh1 = period(c0+xr,c1+yr,c2+zr,NGX,NGY,NGZ)
    mesh2 = period(c0-xr,c1-yr,c2-zr,NGX,NGY,NGZ)
    F[0]=P[mesh1[0]*NGY*NGZ+c1*NGZ+c2,3] - P[mesh2[0]*NGY*NGZ+c1*NGZ+c2,3]
    F[0]=F[0]/distance(P[mesh1[0]*NGY*NGZ+c1*NGZ+c2,:3],P[mesh2[0]*NGY*NGZ+c1*NGZ+c2,:3])
    F[1]=P[c0*NGY*NGZ+mesh1[1]*NGZ+c2,3] - P[c0*NGY*NGZ+mesh2[1]*NGZ+c2,3]
    F[1]=F[1]/distance(P[c0*NGY*NGZ+mesh1[1]*NGZ+c2,:3],P[c0*NGY*NGZ+mesh2[1]*NGZ+c2,:3])
    F[2]=P[c0*NGY*NGZ+c1*NGZ+mesh1[2],3] - P[c0*NGY*NGZ+c1*NGZ+mesh2[2],3]
    F[2]=F[2]/distance(P[c0*NGY*NGZ+c1*NGZ+mesh1[2],:3],P[c0*NGY*NGZ+c1*NGZ+mesh2[2],:3])

    return(F)

def list_2_matrix(Potential,NGX,NGY,NGZ):
    # Convert the linear list of numbers to a grid
    i = 0
    j = 0
    k = 0
    Potential_grid = numpy.zeros(shape=(NGX,NGY,NGZ))
    for i in range (0, NGZ):
     for j in range (0, NGY):
      for k in range (0, NGX):
       Potential_grid[k,j,i] = Potential[k+(j)*NGX+(i)*NGY*NGX]

    return Potential_grid

def grid_2_xyz(P,lattice):
    """Convert the cell shaped grid to pure cartesian"""
    P_xyz = np.zeros(shape=(P.shape[0]*P.shape[1]*P.shape[2],4))
    fa=sqrt(sum(lattice[0,:]**2))/P.shape[0]
    fb=sqrt(sum(lattice[1,:]**2))/P.shape[1]
    fc=sqrt(sum(lattice[2,:]**2))/P.shape[2]

    for i in range (P.shape[0]):
        for j in range (P.shape[1]):
            for k in range (P.shape[2]):
                index = k+j*P.shape[2]+i*P.shape[2]*P.shape[1]
 		P_xyz[index,0] = float(i*fa)
 		P_xyz[index,1] = float(j*fb)
 		P_xyz[index,2] = float(k*fc)
		P_xyz[index,3] = float(P[i,j,k])

    return P_xyz

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
 if i == 2:
  scale_factor = float(inp[0])
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

for i in range(2):
 for j in range(2):
  lattice[i,j] = lattice[i,j]*scale_factor
  
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
 if i > num_atoms + 9 and i < num_atoms + 10 + NGX*NGY*NGZ/5:
  Potential[k]   = inp[0]
  Potential[k+1] = inp[1]
  Potential[k+2] = inp[2]
  Potential[k+3] = inp[3]
  Potential[k+4] = inp[4]
  k = k + 5
  if math.fmod(k,100000) == 0:
   print "Reading potetial, at point", k
  
print ("Average of the potential = ",numpy.average(Potential))
# Start processing the data
average_type=raw_input("Which kind of average would you like? (P)lanar/(S)pherical/(Po)int spherical average ")

#3 Convert the data to a more friendly format

if average_type != "P":
 Potential_grid = list_2_matrix(Potential,NGX,NGY,NGZ)
 Potential_grid = grid_2_xyz(Potential_grid,lattice)
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
 Vacuum_potential = str(Macro_Potential[0])
 Centre_potential = str(Macro_Potential[NGZ/2])
 f.write("Average bulk potential at centre of slab: ")
 f.write(Centre_potential)
 f.write('\n')
 f.write("Average potential of vacuum : vacuum_potential: ")
 f.write(Vacuum_potential)
 f.write('\n')
elif average_type == 'S':
 f = open('Spherical_Av.dat','w')
 f.write('Spherical Average of Potential \n')
 spherical_average=numpy.zeros(shape=(3))
 spherical_av_potential=numpy.zeros(shape=(NGZ,4))
 centroid = numpy.zeros(shape=(3))
 centroid[0] = raw_input("a value on the plane to centre the sphere (Cartesian) ")
 centroid[1] = raw_input("b value on the plane to centre the sphere (Cartesian) ")
 radius = float(raw_input("What radius would you like for spherical averaging? "))
 for i in range (NGZ):
  latt = lattice_vectors(lattice)
  centroid[2] = i*latt[2]/NGZ
  spherical_average = point_sphere(NGX,NGY,NGZ,Potential_grid,axis,centroid,lattice,f,radius)
  print(spherical_average)
  spherical_av_potential[i,0] = i*latt[2]/NGZ
  spherical_av_potential[i,1] = spherical_average[1]
  spherical_av_potential[i,2] = spherical_average[2]
#  spherical_av_potential[i,3] = spherical_average[2]
#  np.savetxt(f,spher_av)
 #print spherical_av_potential
elif average_type == 'Po':
 f = open('Spherical_Av.dat','w')
 f.write('Spherical Average of Potential \n')
 number_points = raw_input("How many points would you like to consider? ")
 for i in range (0,int(number_points)):
  spherical_average=numpy.zeros(shape=(3))
  spherical_av_potential=numpy.zeros(shape=(NGZ))
  centroid = numpy.zeros(shape=(3))
  centroid[0] = raw_input("a lattice value to centre the sphere (Cartesian) ")
  centroid[1] = raw_input("b lattice value to centre the sphere (Cartesian) ")
  centroid[2] = raw_input("c lattice value to centre the sphere (Cartesian) ")
  f.write('Spherical Average of Potential at point \n')
  np.savetxt(f,centroid)
  f.write('\n')
  radius = float(raw_input("What radius would you like for spherical averaging? "))
  spherical_average = point_sphere(NGX,NGY,NGZ,Potential_grid,axis,centroid,lattice,f,radius)
  central_potential = cent_potential(NGX,NGY,NGZ,Potential_grid,centroid,lattice)
  f.write("Centre   Average    Variance \n")
  np.savetxt(f,spherical_average)
  print(" Centre   Average    Variance")
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
 plt.rc('axes', color_cycle=['FF8000', 'CCCC00', '3399FF', '99004C'])


 plt.xlabel('$z / \AA$ ',fontsize=22)
# plt.grid(True)
 if average_type == 'P':
  plt.ylabel('$V_{planar}(z)$ /  $eV$ $\AA^{-1} $',fontsize=22)
  #plt.rc('lines', linewidth=1.5)
  plt.plot(ZAxis,Macro_Potential,'#CC6600')
  plt.plot(ZAxis, Plane_Potential_New,'#0000FF')
  fig = plt.gcf()
  fig.set_size_inches(5,5)
  plt.savefig('Planar.eps',dpi=300)
  plt.savefig('Planar.tif',dpi=450)
  plt.show()
  print("Central bulk potential: ",Centre_potential)
  print("Vacuum potential: ",Vacuum_potential)
 elif average_type == 'S':
  yerr=spherical_av_potential[:,2]
  plt.axhline(ls = '--',linewidth = 1.5, c = 'black')
  #plt.text(max(spherical_av_potential[:,0]/2),'Cell average potential')
  plt.ylabel('$\Phi_{av}(z)$ / $eV \AA^{-3} $',fontsize=22)
  plt.xlim((min(spherical_av_potential[:,0]),max(spherical_av_potential[:,0])))
  plt.errorbar(spherical_av_potential[:,0],spherical_av_potential[:,1],yerr=spherical_av_potential[:,2],lw=1.5,elinewidth=0.5,capsize=3)
  plt.savefig('Error.eps')
  plt.show()
  plt.plot(spherical_av_potential[:,0],spherical_av_potential[:,1],lw=1.5)
  plt.savefig('Normal.eps')
  plt.show()

