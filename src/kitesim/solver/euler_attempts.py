

#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

def explicit_euler(dt,t_end,x0,function,arguments):
    
    # Define time steps
    t = np.arange(0, t_end + dt, dt) # Numerical grid
    # Define the solution space as big as number of time-steps
    x = np.zeros((len(t),len(x0)))
    # Set the initial condition
    x[0] = x0

    # Explicit Euler Method
    for i in range(0, len(t) - 1): # Loop over time-steps
        # update the next solution, using the previous and a function evaluation
        x[i + 1] = x[i] + dt*function(t[i], x[i],arguments)

    return t,x

def implicit_euler(dt, t_end, x0, function,arguments):
    # Define time steps
    t = np.arange(0, t_end + dt, dt)  # Numerical grid
    # Define the solution space as big as the number of time-steps
    x = np.zeros((len(t),len(x0)))
    # Set the initial condition
    x[0] = x0

    for i in range(0, len(t) - 1):  # Loop over time-steps
        # Define the function whose root we want to find
        def root_function(x_next,arguments):
            return x[i] + dt * function(t[i + 1], x_next,arguments) - x_next

        # Use SciPy's root-finding algorithm to find the root
        result = root(root_function, x[i],args=arguments)
        x[i + 1] = result.x[0]

    return t, x

# Define parameters
dt = 0.1 # Time step
t_end = 5 # End time
x0 = np.array([10,0]) # Initial condition
rest_length = x0[0]-x0[1]
spring_constant = 10
damping_constant = 0.1
force_external = -0.5
mass = 10
args =[rest_length, spring_constant, damping_constant, force_external,mass]

def function(t, x, args):
    # the time is fictitious, just there to define a space over which to iterate
    # unpacking arguments
    rest_length, spring_constant, damping_constant, force_external,mass = args
    force = np.zeros(x.shape)
    force[0] = 0
    force[1] = ((x[0]-x[1])-rest_length)*-spring_constant + force_external 
    return force

# Compute solutions
t, explicit_result = explicit_euler(dt, t_end, x0, function,args)
t, implicit_result = implicit_euler(dt, t_end, x0, function,args)

# Plot the results
plt.figure(figsize=(12, 8))
plt.plot(t, explicit_result[:,1], 'bo--', label='Explicit Approximate')
plt.plot(t, implicit_result[:,1], 'ro--', label='Implicit Approximate')
# plt.plot(t, exact_result, 'g', label='Exact')
plt.title('Checking explicit and implicit Euler methods')
plt.xlabel('t')
plt.ylabel('height of lower node')
plt.grid()
plt.legend(loc='lower right')
plt.show()


#%%

import sympy as sp
from sympy import lambdify

class MassSpringDamperNode:
    def __init__(self, m, k, c, F):
        self.m = m  # Mass
        self.k = k  # Spring constant
        self.c = c  # Damping coefficient
        self.F = F  # External force

        self.t = sp.symbols('t')  # Time
        self.x = sp.Function(f'x')(self.t)  # Displacement
        self.v = sp.Function(f'v')(self.t)  # Velocity

        # Define the derivatives
        self.x_dot = self.x.diff(self.t)
        self.v_dot = self.v.diff(self.t)

def create_mass_spring_damper_system(nodes):
    num_nodes = len(nodes)
    variables = []

    for i, node in enumerate(nodes):
        variables.extend([node.m, node.k, node.c, node.F])

    # Define equations of motion
    equations = []

    for i in range(num_nodes):
        m_i, k_i, c_i, F_i = variables[i * 4:i * 4 + 4]
        x_i = nodes[i].x
        v_i = nodes[i].v
        if i == 0:
            # Fixed boundary condition at node 0
            equations.append(m_i * x_i.diff(node.t, 2) - k_i * x_i - c_i * v_i - F_i)
        else:
            x_prev = nodes[i - 1].x
            v_prev = nodes[i - 1].v
            equations.append(m_i * x_i.diff(node.t, 2) - k_i * (x_i - x_prev) - c_i * (v_i - v_prev) - F_i)

    return equations

# Example usage with 3 nodes
node1 = MassSpringDamperNode(m=1.0, k=2.0, c=0.1, F=0.0)
node2 = MassSpringDamperNode(m=2.0, k=1.5, c=0.2, F=0.0)
node3 = MassSpringDamperNode(m=1.5, k=3.0, c=0.15, F=0.0)

nodes = [node1, node2, node3]

# Create the symbolic system with physical properties
equations = create_mass_spring_damper_system(nodes)

# Display the symbolic solution for each node
for i, node in enumerate(nodes):
    print(f"Node {i}:\n{equations[i]}\n")

# Define a function to compute the derivatives of the state variables
def compute_derivatives(state, t, equations):
    derivatives = []
    for i, eqn in enumerate(equations):
        variables = {f'x{i}': state[i] for i in range(len(state))}
        derivatives.append(eqn.subs(variables))
    return derivatives

# Define initial conditions
num_nodes = 3  # Adjust as needed
initial_state = [0.0] * (2 * num_nodes)  # [x1, v1, x2, v2, ..., xn, vn]

# Define time parameters
t_start = 0
t_end = 10
dt = 0.01
timesteps = int((t_end - t_start) / dt) + 1
t = np.linspace(t_start, t_end, timesteps)

# Create lambdified functions for displacement and velocity
displacement_functions = [lambdify(node.t, node.x, "numpy") for node in nodes]
velocity_functions = [lambdify(node.t, node.v, "numpy") for node in nodes]

# Initialize arrays to store the results
displacement_results = [np.zeros(timesteps) for _ in range(num_nodes)]
velocity_results = [np.zeros(timesteps) for _ in range(num_nodes)]

# Simulate the system using Explicit Euler
for i in range(num_nodes):
    state = initial_state.copy()
    for j in range(timesteps):
        # Evaluate displacement and velocity using numerical functions
        displacement_results[i][j] = displacement_functions[i](t[j], state[i * 2])
        velocity_results[i][j] = velocity_functions[i](t[j], state[i * 2 + 1])

        derivatives = compute_derivatives(state, t[j], equations)
        state = [state[k] + dt * derivative for k, derivative in enumerate(derivatives)]

# Plot the results
for i in range(num_nodes):
    plt.plot(t, displacement_results[i], label=f'Node {i + 1} (Explicit Euler)')
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.title('Mass-Spring-Damper System Simulation (Explicit Euler)')
plt.legend()
plt.grid(True)
plt.show()
plt.show()