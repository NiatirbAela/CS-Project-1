import math                                  #trig functions for cosine and sine and radians converter

class LocalDisplacementSolver:               #class 
    def __init__(self, E, A, L, theta_deg, u1, v1, u2, v2):   #stores the inputs
        self.E=E                           #Young's modulus in Pa or Psi
        self.A=A                           #cross-sectional area
        self.L=L                           #element length
        self.theta=theta_deg               #element angle in degrees
        self.u1=u1                         #global x-displacement of node 1
        self.v1=v1                         #global y-displacement of node 1
        self.u2=u2                         #global x-displacement of node 2
        self.v2=v2                         #global y-displacement of node 2

    def compute_local_displacements(self):  #for converting global to local displacements so u to u'
        th=math.radians(self.theta)        #convert degrees to radians
        c, s=math.cos(th), math.sin(th)    #for direction
        self.u1p=c*self.u1+s*self.v1               #local displacement at node 1 which is calculated by u1'=c*u1+s*v1
        self.u2p = c * self.u2 + s * self.v2         #local disp at node 2 
        return self.u1p, self.u2p            #return local displacements

    def compute_strain_stress_force(self):    #to calculate stress, strain and force
        if self.L <= 0:
            raise ValueError("L must be positive.")   #error if number is not positive
        eps=(self.u2p-self.u1p)/self.L         #strain=Δu'/L
        sigma = self.E * eps             #stress=E*strain
        N = sigma * self.A                   #internal force=σA
        return eps, sigma, N
