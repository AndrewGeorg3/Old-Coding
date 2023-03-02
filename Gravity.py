"""
Computer Modelling Ex3: time integration of
n particles moving in a Gravitational potential
using the Velocity Verlet method.

Produces a plot of the total energy of the system
as a function of time. It also calculates the apo/periapses
and the orbital period and saves these to files

The potential is the gravitational potential: V(x) = -GMm/r
where the inputs are taken from files and the objects
are passed through the Particle 3D class
"""

import sys
import math
import numpy as np                 #import the necessary python add ons
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D         #import the particle3D class

#Calculate force and potential energy

def force_Grav(p1, p2, r12):
    """
    Method to return the force on a particle
    in a gravitational potential.
    Force is given by
    F(x) = GMm/r^2

    :param particle: Particle3D instance
    :return: force acting on particle as a float
    """
    mass1 = p1.mass
    mass2 = p2.mass
    G = 1.48814e-34
    force1 = r12*(G*mass1*mass2)/(np.linalg.norm(r12))**3              #Force for a 2 particle system in a Gravitational potential

    return force1

def pot_energy_Grav(p1, p2, r12):
    """
    Method to return potential energy 
    of 2 particles in a gravitational potential
    V(x) = -GMm/r

    :param particle: Particle3D instance
    :return: potential energy of particle as float
    """
    mass1 = p1.mass
    mass2 = p2.mass
    G = 1.48814e-34

    potential = (-1*G*mass1*mass2)/np.linalg.norm(r12)

    return potential            #Potential energy for a 2 particle system in a gravitational potential

def Forces(particlelist):                       #Forces function for looping over n particles 
    forceplanets = []
    for i in range(len(particlelist)):
        netforce = np.array([0.0, 0.0, 0.0])
        for j in range(len(particlelist)):
            if i != j and i > j :                          #doesn't count interactions with itself and gets rid of double copies Fji = Fij
                netforce += force_Grav(particlelist[i], particlelist[j], Particle3D.vector_particle_separation(particlelist[i], particlelist[j]))
        forceplanets.append(netforce)
        
    return np.array(forceplanets)

def Pcorrection():                          #Function accounting for the centre of mass correction
    total_momentum = 0.0                    #this is to stop the system drifting over long run times
    total_mass = 0.0
    
    for i in range(len(particlelist)):
        mass = particlelist[i].mass
        momentum = mass*(particlelist[i].vel)
        total_momentum += momentum
        total_mass += mass
    vcom = total_momentum/total_mass                   #Centre of mass correction to be subtracted from initial velocity
    
    for j in range(len(particlelist)):
        particlelist[j].vel -= vcom                   #changes the initial velocities of all objects
            
def Energy(particlelist):                       #Energy function to loop the gravitaional gravity over all particles
    energy_planets = []
    for i in range(len(particlelist)):
        netenergy = np.array([0.0])
        for j in range(len(particlelist)):
            if i != j:                        #This ignores the case of calculating the potential with itself
                netenergy += 0.5*(pot_energy_Grav(particlelist[i], particlelist[j], Particle3D.vector_particle_separation(particlelist[i], particlelist[j])))
        energy_planets.append(netenergy)             #netenergy is divided by 2 to remove the double counting of interactions

    return np.array(energy_planets)

# Begin main code
def main():

    if len(sys.argv)!=8:          #requires 8 arguments as inputs, will specify which are missing
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <parameters(num and time step)>, <objects file(planets and details)>, <VMD file output (.xyz)>, <Energy output>, <Apoapsis>, <periapsis>, <Orbital periods>")  
        quit()
    else:                             #if the required number of arguments are given, assigns them to the values below
        infile_name   = sys.argv[1]
        object_name   = sys.argv[2]
        xyz_filename  = sys.argv[3] 
        energy_filename = sys.argv[4]
        largepositions_filename = sys.argv[5]
        smallpositions_filename = sys.argv[6]
        orbitalperiods_filename = sys.argv[7]
        
    input_handle = open(infile_name, "r") 
    input_parameters = input_handle.readline().split(',')          #Splits the parameters file with a comma
    dt = float(input_parameters[0])                                #assigning the parameters data to our variables
    numstep = int(input_parameters[1])
    input_handle.close()

    # Open output files
    outfile_traj = open(xyz_filename, "w")
    outfile_Verlet_E = open(energy_filename, "w")
    outfile_Apses = open(largepositions_filename, "w")
    outfile_periapses = open(smallpositions_filename, "w")
    outfile_orbits = open(orbitalperiods_filename, "w")

    # Set up simulation parameters
    time_Verlet = 0.0

    readfile = open(object_name, "r")
    particlelist = Particle3D.read_from_file(readfile)      #n particles are now taken from my Particle3D

             
    #write out initial conditions for Energy
    for i in range(len(particlelist)):
        energy_Verlet = np.linalg.norm(particlelist[i].kinetic_energy() + Energy(particlelist))

    outfile_Verlet_E.write("{0:f} {1:12.8f}\n".format(time_Verlet, energy_Verlet))

    # Initialise data lists for plotting later
    time_list_Verlet = [time_Verlet]
    energy_list_Verlet = [energy_Verlet]
        
    orbit_list = []                       #Attempted method for calculating the orbital period- Unsuccessful
    for i in range(len(particlelist)):
        xvel = np.linalg.norm(particlelist[i].velocity)
        orbit_list.append(xvel)	
        
    max_list = []                        #Method for calculating the apoapses.
    for i in range(len(particlelist)):
        max_list.append(0.0)                 #starts at 0 as the separation can only be larger
    
    min_list = []                               #Method for calculating the periapses
    for i in range(len(particlelist)):          #starts at the initial separation as the periapsis can only be smaller than this
        min_list.append(np.linalg.norm(Particle3D.vector_particle_separation(particlelist[i], particlelist[0])))

    # Start the time integration loop for Velocity Verlet
    for n in range(numstep):
    
        forceplanets = Forces(particlelist)
    
        for j in range(len(particlelist)):
            
            # Update particle position
            particlelist[j].leap_pos2nd(dt, forceplanets[j])

        force_new = Forces(particlelist)

            # Update particle velocity by averaging
            # current and new forces
        for k in range(len(particlelist)):
            
            particlelist[k].leap_velocity(dt, 0.5*(forceplanets[k] + force_new[k]))
            new_velocity = [particlelist[k].velocity]
            
            xvel = np.linalg.norm(new_velocity[0])

            if (xvel > orbit_list[k]):
                orbit_list[k] = time_Verlet

        """
        Above was an attempt made for orbital period.
        The logic was that by taking the velocity in 1 dimension, 
        the max value would occur every time 1 orbit was completed.
        By saving the numberstep that this occured at
        (for every planet being looped over), the orbital 
        period would be able to be saved and displayed in an outfile.
        This seems to work for some examples, but is not consistent
        and we could not work out why.
        Other methods were explored, such as the turning point 
        (value going from positive to negative) and using radians
        for a percentage of an orbit. However, these were 
        unable to be coded either
        """

            # Re-define force value
        for l in range(len(particlelist)):
             forceplanets[l] = force_new[l]

        # Increase time
        time_Verlet += dt

        #update particle separation magnitude
        for m in range(len(particlelist)):
            for p in range(len(particlelist)):
                if m != p:                        #doesn't count interactions with itself
                    r12 = Particle3D.vector_particle_separation(particlelist[p], particlelist[m])
                    R12 = np.linalg.norm(r12)
                    
        #loop the total energy
        for z in range(len(particlelist)):
            energy_Verlet = np.linalg.norm(Energy(particlelist) + particlelist[z].kinetic_energy())
            
        for i in range(len(particlelist)):
            if i != 0:                                          #if the particle is not the 0th element (the sun), calculate the separation
                R12 = np.linalg.norm(Particle3D.vector_particle_separation(particlelist[i], particlelist[0]))        #only count separations beteween bodies and the sun
                if (particlelist[i].label == "Moon"):                     #special case for the moon- only consider the orbit around the earth
                    R12 = np.linalg.norm(Particle3D.vector_particle_separation(particlelist[i], particlelist[3]))

                if (R12 > max_list[i]):               #everytime R12 is greater than the previous value, update the stored value in the list 
                    max_list[i] = R12                 # at the end of the loop, this is then saved to the outfile for the Apses
                    
                if (R12 < min_list[i]):               #same process for periapses
                    min_list[i] = R12
                    
        
        # Output particle information and write to file

        if n%10==0:            #selects every 10th step to vmd outfile
            outfile_traj.write("{0:n}\n".format(len(particlelist)))
            outfile_traj.write(f"point = {time_Verlet}\n") #Joe says this is a cool new way of formatting python strings
            for i in range(len(particlelist)):
                outfile_traj.write(str(particlelist[i]))
                outfile_traj.write("\n")


        outfile_Verlet_E.write("{0:f} {1:12.8f}\n".format(time_Verlet, energy_Verlet))          #write the energy to an outfile

        # Append information to data lists
        time_list_Verlet.append(time_Verlet)
        energy_list_Verlet.append(energy_Verlet)

    for i in range(len(particlelist)):
        outfile_Apses.write("{0:s} {1:12.8f}\n".format(particlelist[i].label, max_list[i]))            #writing the outputs to a file, 
        outfile_periapses.write("{0:s} {1:12.8f}\n".format(particlelist[i].label, min_list[i]))        #by looping over all particles 
        outfile_orbits.write("{0:s} {1:12.8f}\n".format(particlelist[i].label, orbit_list[i]))
        
    # Post-simulation:
    # Close output file
    outfile_Verlet_E.close()

    #Plot particle energy to screen (Verlet)
    pyplot.title('Velocity Verlet: Total Energy vs Time')
    pyplot.xlabel('Time (Days)')
    pyplot.ylabel('Energy (kg Au Day^(-2)')
    pyplot.plot(time_list_Verlet, energy_list_Verlet, label = "Total Energy")
    pyplot.legend(loc='upper left')
    pyplot.show()

 # Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
