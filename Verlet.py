"""
Computer Modelling Ex2: time integration of
two particles moving in a Morse potential
using Symplectic Euler and Velocity Verlet methods.
This code details Velocity Verlet

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to files

The potential is the Morse potential: V(x) = D((1- exp(-a(r12-re)))(^2)-1)
where a, D, r2, r1 and re are constants or inputs taken from my
Particle 3D class
"""

import sys
import math
import numpy as np                 #import the necessary python add ons
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D         #import the particle3D class

#Calculate force and potential energy

def force_Morse(particle1, particle2, a, r12, R12, re, D):
    """
    Method to return the force on a particle
    in a Morse potential.
    Force is given by
    F(x) = -dV/dx = 2*a*D*(1-exp(-a(R12-re)))*exp(-a(R12-re))*((r12)/(R12))

    :param particle: Particle3D instance
    :return: force acting on particle as a float
    """
    r12 = Particle3D.vector_particle_separation(particle1, particle2)     #particle separation
    R12  = np.linalg.norm(r12)                           #magnitude of the particle separation
    unitSep = r12/R12
    E = (math.exp(-a*(R12-re)))                         #the exponential term for the force
    force1 = 2*a*D*(1-E)*(E)*unitSep              #Force for a 2 particle system in a Morse potential
    return force1


def pot_energy_Morse(particle1, particle2, R12, D, re, a):
    """
    Method to return potential energy 
    of 2 particles in a Morse potential
    V(x) = D((1-exp(-a*(R12-re)))(^2)-1)

    :param particle: Particle3D instance
    :return: potential energy of particle as float
    """
    r12 = Particle3D.vector_particle_separation(particle1, particle2)
    R12  = np.linalg.norm(r12)
    E = (math.exp(-a*(R12-re)))
    potential = D*(((1-E)**2)-1)
    return potential            #Potential energy for a 2 particle system in a Morse potential


# Begin main code
def main():

    """
    file_handle = open("initial_values_test",'r')
    p2 = Particle3D.read_from_file(file_handle)
    p1 = Particle3D ("particle1", -0.5, 1.0, 0.0, 0.1, 0.0, 0.0, 2.0)
    print(p1)
    print(p2)
    print(p1.position)
    print(p2.position)
    print(Particle3D.vector_particle_separation(p1,p2))
    print(p2.position-p1.position)
    # Read name of output file from command line
    """

    if len(sys.argv)!=6:          #requires 6 arguments as inputs, will specify which are missing
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file>, <time step>, <num step>, <input file(Oxygen OR Nitrogen>, <constants file(Oxygen OR Nitrogen constants>")  
        quit()
    else:                             #if the required number of arguments are given, assigns them to the values below
        outfile_name = sys.argv[1]              #file output name
        dt = float(sys.argv[2])                 #time step
        numstep = int(sys.argv[3])              #number step
        infile_name = sys.argv[4]               #Particle to be simulated- either Oxygen or Nitrogen
        constants_name = sys.argv[5]            #the required constants for the chosen particle

    # Open output file
    outfile_Verlet = open(outfile_name, "w")

    # Set up simulation parameters
    time_Verlet = 0.0

    readfile = open(infile_name, "r")             #reads a particle input file as named in the terminal
    constantsfile = open(constants_name, "r")           #reads the file for constants specified in the terminal
    line = constantsfile.readline().split(',')           #splits the constants by a comma
    D = float(line[0])                            #assign the values from the data file to my constants
    re = float(line[1])
    a = float(line[2])

    p1, p2 = Particle3D.read_from_file(readfile)        #particle 1 and 2 are now taken from my Particle3D

    r12 = Particle3D.vector_particle_separation(p1, p2)    #vector particle separation
    R12  = np.linalg.norm(r12)                             #magnitude of particle separation    

    force1 = force_Morse(p1, p2, a, r12, R12, re, D)     #define force 1 before the loop so there is a start point
    force2 = -force1

    # Write out initial conditions for Verlet

    energy_Verlet = p1.kinetic_energy() + p2.kinetic_energy()+ pot_energy_Morse(p1, p2, R12, D, re, a)

    #writes the initial conditions to my file

    outfile_Verlet.write("{0:f} {1:f} {2:12.8f}\n".format(time_Verlet, R12, energy_Verlet))

    # Initialise data lists for plotting later
    time_list_Verlet = [time_Verlet]
    pos1_list_Verlet = [p1.position]
    pos2_list_Verlet = [p2.position]
    energy_list_Verlet = [energy_Verlet]
    sep_list_Verlet = [R12]

    # Start the time integration loop for Velocity Verlet
    for i in range(numstep):

        # Update particle position
        p1.leap_pos2nd(dt, force1)
        p2.leap_pos2nd(dt, force2)

        # Update force
        force_new = force_Morse(p1, p2, a, r12, R12, re, D)
        force_new2 = -force_new

        # Update particle velocity by averaging
        # current and new forces
        p1.leap_velocity(dt, 0.5*(force1+force_new))
        p2.leap_velocity(dt, 0.5*(force2+force_new2))

        # Re-define force value
        force1 = force_new
        force2 = -force_new

        # Increase time
        time_Verlet += dt

        #update particle separation magnitude
        R12 = np.linalg.norm(Particle3D.vector_particle_separation(p1, p2))

        # Output particle information and write to file
        energy_Verlet = p1.kinetic_energy() + p2.kinetic_energy() + pot_energy_Morse(p1, p2, R12, D, re, a)
        outfile_Verlet.write("{0:f} {1:f} {2:12.8f}\n".format(time_Verlet, R12, energy_Verlet))

        # Append information to data lists
        time_list_Verlet.append(time_Verlet)
        pos1_list_Verlet.append(np.linalg.norm(p1.position))
        pos2_list_Verlet.append(np.linalg.norm(p2.position))
        energy_list_Verlet.append(energy_Verlet)
        sep_list_Verlet.append(R12)

    # Post-simulation:
    # Close output file
    outfile_Verlet.close()

    # Plot particle trajectory to screen (Verlet)
    pyplot.title('Velocity Verlet: Separation vs Time')
    pyplot.xlabel('Time (s x10^-14)')
    pyplot.ylabel('Particle Separation (Å)')
    pyplot.plot(time_list_Verlet, sep_list_Verlet)
    pyplot.show()

    # Plot particle energy to screen (Verlet)
    pyplot.title('Velocity Verlet: Total Energy vs Time')
    pyplot.xlabel('Time (s x10^-14)')
    pyplot.ylabel('Energy (eV)')
    pyplot.plot(time_list_Verlet, energy_list_Verlet)
    pyplot.show()

# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
