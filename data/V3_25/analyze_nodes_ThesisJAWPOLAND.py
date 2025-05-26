import csv
import math

# Path to the CSV file
csv_file = "nodes_ThesisJAWPOLAND.csv"

# Indices to compute distances between (as tuples)
index_pairs = [
    (0, 26),
    (0, 21),
    (21, 27),
    (26, 27),
    (26, 20),
    (27, 28),
    (27, 29),
    (29, 14),
    (29, 13),
    (28, 12),
    (28, 11),
    (0, 33),
    (33, 34),
    (33, 20),
    (27, 20),
    (34, 9),
    (34, 8),
    (33, 35),
    (35, 6),
    (35, 7),
]

# Read the CSV and extract relevant columns
points = []
with open(csv_file, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Extract coordinates as floats
        x = float(row["x-coordinate [mm]"])
        y = float(row["y-coordinate [mm]"])
        z = float(row["z-coordinate [mm]"])
        points.append((x, y, z))


# Function to compute Euclidean distance
def distance(p1, p2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


# Print distances
for i, j in index_pairs:
    d = distance(points[i], points[j])
    print(f"Distance between indices {i} and {j}: {d:.3f} mm")
