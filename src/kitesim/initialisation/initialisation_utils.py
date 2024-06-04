import sys
import importlib.util
from scipy.spatial import ConvexHull
import numpy as np


def load_module_from_path(module_name, module_path):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def calculate_projected_area(points):
    # Project points onto the x,y plane
    xy_points = points[:, :2]

    # Find the convex hull
    hull = ConvexHull(xy_points)
    hull_points = xy_points[hull.vertices]

    # Using the shoelace formula
    x = hull_points[:, 0]
    y = hull_points[:, 1]

    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
