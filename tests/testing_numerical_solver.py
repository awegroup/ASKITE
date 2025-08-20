# # %%
# import sympy as sp

# # Define the symbols
# A, B, C, D, E, x = sp.symbols("A B C D E x")

# # Define specific values for A, B, C, D, and E
# A_value = 22  # m
# B_value = 28  # v
# C_value = 2000  # Fa,x
# D_value = 200  # Lt
# E_value = 8000  # Fa,z

# # Define the equation
# equation = sp.Eq(A * (B**2 / x), C * sp.tan(sp.asin(x / D)) + E)


# # Substitute the specific values
# equation_with_values = equation.subs(
#     {A: A_value, B: B_value, C: C_value, D: D_value, E: E_value}
# )

# # Solve the equation for x
# solutions = sp.solve(equation_with_values, x)

# # # Print the solutions
# # print(solutions)

# # Substitute the specific values into the solutions and evaluate them numerically
# numerical_solutions = [
#     sol.evalf(subs={A: A_value, B: B_value, C: C_value, D: D_value, E: E_value})
#     for sol in solutions
# ]

# # Print the numerical solutions
# print(numerical_solutions)

# for sol in numerical_solutions:
#     if sp.re(sol) >= 0:
#         solution_sympy = sp.re(sol)
#         print(sp.re(sol))

# # %% [markdown]
# # We now have a solution, let's try to Wolfram Alpha's solution

# # %%
# # Define the equation
# expression = (A * E * B**2) / (2 * (C**2 + E**2)) - (1 / 2) * sp.sqrt(
#     (A**2 * E**2 * B**4) / (C**2 + E**2) ** 2
#     - (2 * (A**2 * B**4 - D**2 * E**2)) / (3 * (C**2 + E**2))
#     + (
#         A**6 * B**12
#         + 36 * A**4 * C**2 * D**2 * B**8
#         - 3 * A**4 * D**2 * E**2 * B**8
#         + 3 * A**2 * D**4 * E**4 * B**4
#         + 18 * A**2 * C**2 * D**4 * E**2 * B**4
#         - D**6 * E**6
#         + 6
#         * sp.sqrt(3)
#         * sp.sqrt(
#             A**10 * C**2 * D**2 * B**20
#             + 8 * A**8 * C**4 * D**4 * B**16
#             - 3 * A**8 * C**2 * D**4 * E**2 * B**16
#             + 16 * A**6 * C**6 * D**6 * B**12
#             + 3 * A**6 * C**2 * D**6 * E**4 * B**12
#             + 20 * A**6 * C**4 * D**6 * E**2 * B**12
#             - A**4 * C**2 * D**8 * E**6 * B**8
#             - A**4 * C**4 * D**8 * E**4 * B**8
#         )
#     )
#     ** (1 / 3)
#     / (3 * (C**6 + 3 * E**2 * C**4 + 3 * E**4 * C**2 + E**6) ** (1 / 3))
#     + (
#         (
#             A**4 * B**8
#             - 12 * A**2 * C**2 * D**2 * B**4
#             - 2 * A**2 * D**2 * E**2 * B**4
#             + D**4 * E**4
#         )
#         * (C**6 + 3 * E**2 * C**4 + 3 * E**4 * C**2 + E**6) ** (1 / 3)
#     )
#     / (
#         3
#         * (C**2 + E**2) ** 2
#         * (
#             A**6 * B**12
#             + 36 * A**4 * C**2 * D**2 * B**8
#             - 3 * A**4 * D**2 * E**2 * B**8
#             + 3 * A**2 * D**4 * E**4 * B**4
#             + 18 * A**2 * C**2 * D**4 * E**2 * B**4
#             - D**6 * E**6
#             + 6
#             * sp.sqrt(3)
#             * sp.sqrt(
#                 A**10 * C**2 * D**2 * B**20
#                 + 8 * A**8 * C**4 * D**4 * B**16
#                 - 3 * A**8 * C**2 * D**4 * E**2 * B**16
#                 + 16 * A**6 * C**6 * D**6 * B**12
#                 + 3 * A**6 * C**2 * D**6 * E**4 * B**12
#                 + 20 * A**6 * C**4 * D**6 * E**2 * B**12
#                 - A**4 * C**2 * D**8 * E**6 * B**8
#                 - A**4 * C**4 * D**8 * E**4 * B**8
#             )
#         )
#         ** (1 / 3)
#     )
# )

# # Substitute the specific numerical values into the equation and evaluate it numerically
# numerical_result = expression.subs(
#     {A: A_value, B: B_value, C: C_value, D: D_value, E: E_value}
# ).evalf()

# # Print the numerical result
# print(numerical_result)

# # %% [markdown]
# # Let's try to solve the equation directly like Wolfram's alphas

# # %%
# import numpy as np
# from scipy.optimize import root


# # Define the function to find roots for
# def equation_to_solve(
#     r_0, m_kite, v_kite, length_tether, force_aero_xwind, force_aero_zwind
# ):
#     return (
#         m_kite * (v_kite**2 / r_0)
#         - force_aero_xwind * np.tan(np.arcsin(r_0 / length_tether))
#         - force_aero_zwind
#     )


# # Defining the variables
# m_kite = 22  # kg
# v_kite = 28  # m/s
# length_tether = 200  # m
# force_aero_xwind = 8000  # N
# force_aero_zwind = 3000  # N

# # Initial guess for the root
# r_0_initial_guess = 20

# # Find the root numerically and get full output
# result = root(
#     equation_to_solve,
#     r_0_initial_guess,
#     args=(m_kite, v_kite, length_tether, force_aero_xwind, force_aero_zwind),
#     method="hybr",
# )

# # Print the full output
# print(f"nfev:{result.nfev}")
# print("Numerical Solution:", result.x)
# print(f"Ratio Fx/Fz = {force_aero_xwind/force_aero_zwind:.2f}")


# # %%
