
#%% 
import numpy as np

## Algorithm (adapted from P.Thedens)

##Setting parameters
""" from Paul's thesis: It should be noted that this is a lumped Massn matrix where 
all off-diagonal terms are zero, and this format allows for a cheap inversion."""
Massn                   = np.array([[1,0],[0,1]]) # example here: 2 nodes, both with Massn = 1
dt                      = 1e-3 # time step
K_spring                = 1e2 # spring stiffness
boolean_external_model  = False # we do it without external model (aerodynamics) here True if aerodynamic model is used, False if not
rest_length             = .9 #spring is prestretched to 1m

##Initialise for current time step _t
# positions
Xn_ini              = np.array([[0,0,0],[0,0,-1]]) # example here: 2 nodes, 1 at origin, 1 at [0,0,-1] (could be a hanging Massn)
Xn                  = Xn_ini

# velocities
Vn_ini              =  np.zeros((2,3))
nodes_vel_t         = Vn_ini
Vn_min_half         = Vn_ini # initialise velocity at t_min_half = t_ini
# kinetic energy Wkin 
Wkin_n              = 1
Wint_n              = 1
# internal forces
F_internal_ini      =  np.zeros((2,3))   #zero internal forces at t = 0
F_internal          = F_internal_ini #initialise internal forces
# external forces
F_external_ini      = np.zeros((2,3)) 
F_external_ini[1,2] = -9.86*1#gravity force at t = 0
F_external          = F_external_ini #gravity force at t = 0
# residual forces
F_residualn          = F_external-F_internal


''' REACTION FORCE
Bit of explanation on reaction forces, and why they are used in a system with constraints:
(taken from: https://manuals.dianafea.com/d103/Theory/Theoryse453.html)

To check the quality of the numerical solution the residual forces can be output. 
The calculated residual forces should be zero or small compared to the external 
and internal forces in case of equilibrium for unconstrained degrees of freedom.

For constrained degrees of freedom (supports) the residual force will have a finite value 
equal to the force that the element model exhibits on the constrained degree of freedom. 
The reaction forces are defined as the forces that constrained degrees of freedom (supports) 
exhibits on the element model. Hence, the reaction forces are opposite to the residual forces 
in the constrained degrees of freedom:
F_reaction = -F_residualn '''
F_reaction      = np.zeros((2,3)) #initialise reaction forces
F_reaction[0]   = 1e-3*np.ones((1,3)) #initialise this as way to low

# the energy tolerance is a ratio of kinetic to internal energy
# which should become lower than the set tolerance, such that all energy is converted to internal energy
tol_energy = 1e-5
# the force tolerance is the difference between the external and internal forces
# where for a convergence system, the residual forces are exactly opposing the reaction forces
tol_force = 1e-5

# print (f"Initial position: {Xn_ini}")
# print (f"Initial velocity: {Vn_ini}")
# print (f"Initial internal force: {F_internal_ini}")
# print (f"Initial external force: {F_external_ini}")
# print (f"Initial residual force: {F_residualn}")
# print (f"Initial reaction force: {F_reaction}")
# print (f"Initial kinetic energy: {Wkin_n}")
# print (f"Initial internal energy: {Wint_n}")
# print (f"Wkin_n/Wint_n: {Wkin_n/Wint_n}")
# print( f"abs(F_residualn)/abs(F_reaction): {np.sum(np.linalg.norm(F_residualn))/(np.sum(np.linalg.norm(F_reaction)))} \n")


#1 Initiliase positions: Xn, Velocities: Vn_plus_half: Kinetic energy Wkin_n
Xn              = Xn_ini
Vn_plus_half    = Vn_ini
Vn_min_half     = Vn_ini
Wkin_n          = 1

i = -1
imax = 2e5
#2 Main loop, which checks for energy and force convergence
while i<imax and (Wkin_n/Wint_n) > tol_energy and (np.sum(np.linalg.norm(F_residualn))/(np.sum(np.linalg.norm(F_reaction)))) > tol_force:

    #3 Initialise kinetic energy Wkin for previous time steps (t-1/2) #eq 3.74
    Vn_plus_half     = (dt/2)*np.dot(np.linalg.inv(Massn),F_residualn)
    Wkin_n_plus_half = np.sum(np.dot(np.dot(Vn_plus_half.T,Massn),Vn_plus_half)) # eq. 3.71 ##TODO: verify if this is correct
    Xn_plus_1        = Xn + (dt/2) * Vn_plus_half
    Vn               = 0 #set to zero at the energy peak Vn_plus_half*(dt/2) # rewritten eq. 3.68 and combined with eq. 3.74 insights

    #4 Determine external force vector fext (Eq. 3.46, 3.52 and 3.56)
    F_external = F_external_ini #gravity force remains constant here
    
    boolean_KDR = False

    #5 KDR loop, which checks for kinetic energy peaks, and loops untill an energy peak has been reached
    while not boolean_KDR and i<imax:
        i += 1
        # print( f"it:{i}-------")# | tol_E: {Wkin_n/Wint_n:.2f} | tol_F: {np.sum(np.linalg.norm(F_residualn))/(np.sum(np.linalg.norm(F_reaction))):.2f}")

        ## Start a new-integration step: (n+1)-->(n) & (n+1/2)-->(n-1/2)
        Vn_min_half       = Vn_plus_half
        Wkin_n_min_half   = Wkin_n_plus_half
        Xn                = Xn_plus_1
        Vn_min_1          = Vn
        F_residualn_min_1 = F_residualn
        Massn_min_1       = Massn

        #6 Determine internal force vector F_internal (Eq. 3.57 and 3.47)
        elongation =  np.linalg.norm(Xn[1]-Xn[0]) - rest_length
        F_internal[0][2] = K_spring*elongation #use Hook's law here
        F_internal[1][2] = -K_spring*elongation #use Hook's law here

        #6 Determine new (fake) Massn-matrix (Eq. 3.76)
        Massn = Massn ##TODO: make it varying and not equal to the physical gravity Massn

        #7 Compute nodal velocity at t+1/2 (Eq. 3.70)
        Vn_plus_half = Vn_min_half + dt * np.dot(np.linalg.inv(Massn),F_residualn)

        #7 Compute nodal position at t+1 (Eq. 3.69)
        Xn_plus_1 = Xn + dt * Vn_plus_half

        #8 Correct Dirichlet boundary conditions 
        Vn_plus_half[0] = Vn_ini[0] # set-vel to zero (initial), to keep the first node at the origin
        Xn_plus_1[0]    = Xn_ini[0] # set-pos to zero (initial), keep the first node at the origin

        #9 if considering changes in the external forces
        if boolean_external_model:

            #10 Update external force vector fext (Eq. 3.56)
            F_external = F_external_ini #gravity force remains constant here 
        
        #12 Compute kinetic energy (Eq. 3.71)
        Wkin_n_plus_half = np.sum(np.dot(np.dot(Vn_plus_half.T,Massn),Vn_plus_half)) # eq. 3.71 ##TODO: verify if this is correct

        #12 Update from previous time steps (kinetic-energy, internal-energy, F_residualn) 
        Vn = Vn_plus_half*dt + Vn_min_1 # rewritten eq. 3.68
        Wkin_n = np.sum(np.dot(np.dot(Vn.T,Massn),Vn))
        
        elongation =  np.linalg.norm(Xn[1]-Xn[0]) - rest_length #determined at n not at n+1
        Wint_n  = 0.5*K_spring*(elongation)**2 # taken from potential energy theory
        
        F_residualn = F_external-F_internal #this is still R(u_n) (i.e. not a function of n+1 but of n, as it should be)

        # print(f"F_residualn (z-component): {F_residualn[:,2]}")
        print(f'Wkin_n+half: {Wkin_n_plus_half:.3f} | Wkin_n-min: {Wkin_n_min_half:.3f} | Wkin_n: {Wkin_n:.5f} | Vn_plus_half: {Vn_plus_half[1,2]:.3f}')
        
        #13 if energy peak is detected (if kinetic energy is decreasing) 
        if Wkin_n_plus_half < Wkin_n_min_half: # then we get a better estimations of the deformation during this peak
            
            print(f'Energy peak detected at iteration {i}')

            #14 Correct deformation
            #first recompute Fresidual at Xn_min_1 (could also be stored from an earlier run...)
            Vn_min_32 = Vn_min_half +dt*np.dot(np.linalg.inv(Massn_min_1),F_residualn_min_1) #important to take old residual here
            Wkin_n_min_32 = np.sum(np.dot(np.dot(Vn_min_32.T,Massn_min_1),Vn_min_32)) # adopted from eq. 3.71
            q = (Wkin_n_min_half - Wkin_n_plus_half ) / ( 2*Wkin_n_min_half - Wkin_n_min_32 - Wkin_n_plus_half ) #eq 3.73 
            Xn = -(1+q)*Vn_min_half + 0.5*(q*dt)*np.dot(np.linalg.inv(Massn),F_residualn_min_1) #eq. 3.72

            #15 Update internal force vector fint (Eq. 3.57 and 3.47)
            elongation =  np.linalg.norm(Xn[1]-Xn[0]) - rest_length
            F_internal[0][2] = K_spring*elongation #use Hook's law here
            F_internal[1][2] = -K_spring*elongation #use Hook's law here

            #Update Massn-matrix (Eq. 3.76)
            Massn = Massn  ##TODO: make it varying and not equal to the physical gravity Massn

            #16 Correct the Dirichlet boundary conditions
            Xn[0] = Xn_ini[0] #keep the first node at the origin
              
            #17 Determine reaction force vector fraction (Eq. 3.75)
            # print(f"Reaction force: {F_reaction}")
            # print(f"Internal force: {F_internal}") 
            F_reaction[0][2] = F_internal[0][2] #reaction force is opposite to internal force at the Dirichlet boundary condition

            #18 Stop loop as kinetic energy is decreasing
            boolean_KDR = True #stop the KDR loop, i.e. set kinetic-energy to zero
        
        #19 end if statement

    #20 end while-loop (KDR loop)
#21 end while-loop (main loop)

print(f"Done in {i} iterations")
print(f"Wkin_n/Wint_n: {Wkin_n/Wint_n} | tol_E: {tol_energy}")
print(f"F_residual/F_reaction: {np.sum(np.linalg.norm(F_residualn))/(np.sum(np.linalg.norm(F_reaction)))} | tol_F: {tol_force}")
print(f"Initial z-position: {(Xn_ini[:,2])}")
print(f"Final z-position  : {(Xn[:,2])}")

# %%
