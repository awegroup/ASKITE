import numpy as np
import matplotlib.pyplot as plt


def generate_distribution(n_points):
    x = np.linspace(0, 1, n_points)
    y = np.zeros_like(x)
    peak_index = int(len(x) * 0.25)
    y[:peak_index] = x[:peak_index] / (0.25)
    y[peak_index:] = (1 - x[peak_index:]) / (0.75)
    return x, y / np.sum(y)


def print_interpolated_distribution(distribution):
    print("Interpolated Distribution:")
    for i, value in enumerate(distribution):
        print(f"x={i/n_points:.2f}: {value:.2f}")


n_points = 50
x, distribution = generate_distribution(n_points)

plt.plot(x, distribution, marker="o")
plt.title("Distribution")
plt.xlabel("x")
plt.ylabel("Probability Density")
plt.grid(True)
plt.show()

print_interpolated_distribution(distribution)
