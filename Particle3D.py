"""
Particle 3D- stage 1 of assignment 2
Defines the classes that will be required in my
main code and creates 2 particles from an input
file. This file will be asked for in my main codes
and the static methods here will create those particles
"""

import numpy as np                #import numpy

class Particle3D(object):

    """
    Class to describe 3D particles.

    Properties:
    position - position along all axes (x,y,z)
    velocity - velocity along all axes (v1,v2,v3)
    mass(float) - particle mass

    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    """

    def __init__(self, label, pos, vel, mass):

        """
        Initialise a Particle3D instance
        """

        self.position = pos
        self.velocity = vel             #defines the classes for all of my particle components
        self.mass = mass
        self.label = label
    

    def __str__(self):

        """
        Define output format.
        For particle p=((2.0, 1.0, 7.0), (0.5, 0.3, 0.1), 1.0) this will print as
        "position = 2.0, 1.0, 7.0" 
        "velocity = 0.5, 0.3, 0.1"
        "mass = 1.0"
        """
        mystring = "{0} {1:12e} {2:12e} {3:12e}".format(self.label, self.position[0], self.position[1], self.position[2])
        return mystring
       
    
    def kinetic_energy(self):

        """
        Return kinetic energy as
        1/2*mass*((v1)^2 + (v2)^2 + (v3)^2)
        """

        return 0.5*self.mass*(np.linalg.norm(self.velocity))**2
        

    # Time integration methods

    def leap_velocity(self, dt, force):

        """
        First-order velocity update,
        v1(t+dt) = v1(t) + dt*F(t)
        v2(t+dt) = v2(t) + dt*F(t)
        v3(t+dt) = v3(t) + dt*F(t) 

        :param dt: timestep as float
        :param force: force on particle as float
        """

        self.velocity += dt*force/self.mass


    def leap_pos1st(self, dt):

        """
        First-order position update,
        x(t+dt) = x(t) + dt*v1(t) 
        y(t+dt) = y(t) + dt*v2(t)  
        z(t+dt) = z(t) + dt*v3(t)

        :param dt: timestep as float
        """

        self.position += dt*self.velocity            #adds the change in position due to velocity to the ne position of a particle


    def leap_pos2nd(self, dt, force):

        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v1(t) + (1/2*(mass))*dt^2*F(t)/mass
        y(t+dt) = y(t) + dt*v2(t) + (1/2*(mass))*dt^2*F(t)/mass
        z(t+dt) = z(t) + dt*v3(t) + (1/2*(mass))*dt^2*F(t)/mass
        :param dt: timestep as float
        :param force: current force as float
        """

        self.position += dt*self.velocity + 0.5*dt**2*force/self.mass
        
    @staticmethod
    def read_from_file(filehandle):
        
        particlelist = []
        for line in filehandle:
           # filehandle.readline()              #reading line 1 of my file input
        
            data = line.split(',')                   #splits the data up by every comma

            pos = np.array([float(data[1]), float(data[2]), float(data[3])])     #assigns the elements of the inputs to the position in x,y and z
            vel = np.array([float(data[4]), float(data[5]), float(data[6])])     #assigns the elements of the inputs to the velocity in v1,v2 and v3
            mass = float(data[7])                    #assigns the mass to be a float
            label = str(data[0])                     #assigns a label to the particle
            particle = Particle3D(label, pos, vel, mass)         #creates particle 1 using the classes above
            particlelist.append(particle)
        return particlelist

    @staticmethod 
    def vector_particle_separation(p1, p2): 
        return p2.position - p1.position                    #vector separation of the particles is the difference in their positions
