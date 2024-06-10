import sys
import importlib.util
from scipy.spatial import ConvexHull
import numpy as np
import attr


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


def update_dict_with_instantiated_classes(
    data_dict: dict, dict_nested_instances: dict
) -> dict:
    """Update the data_dict to include the instantiated classes.

    Args:
        data_dict (dict): The dictionary with the data to be used to create the class.
        dict_nested_instances (dict): The dictionary with the instantiated classes.

    Returns:
        The updated data_dict.
    """
    for key, value in data_dict.items():
        if key in dict_nested_instances:
            data_dict[key] = dict_nested_instances[key]
    return data_dict


def create_attr_class_from_dict(class_name: str, nest0_data_dict: dict):
    """Create an attrs class from a dictionary with up to three levels of nesting.

    Args:
        class_name (str): The name of the class to be created.
        nest0_data_dict (dict): The dictionary with the data to be used to create the class.

    Returns:
        The instantiated class.
    """
    nest0_attributes = {}
    nest1_attributes = {}
    nest2_attributes = {}
    nest3_attributes = {}
    nest4_attributes = {}

    nest1_instances = {}
    nest2_instances = {}
    nest3_instances = {}
    nest4_instances = {}

    # Level 0
    for nest0_key, nest0_value in nest0_data_dict.items():
        if isinstance(nest0_value, dict):
            nest1_data_dict = nest0_value

            # Level 1
            for nest1_key, nest1_value in nest1_data_dict.items():
                if isinstance(nest1_value, dict):
                    nest2_data_dict = nest1_value

                    # Level 2
                    for nest2_key, nest2_value in nest2_data_dict.items():
                        if isinstance(nest2_value, dict):
                            nest3_data_dict = nest2_value

                            # Level 3
                            for nest3_key, nest3_value in nest3_data_dict.items():
                                if isinstance(nest3_value, dict):
                                    nest4_data_dict = nest3_value

                                    # Level 4
                                    for (
                                        nest4_key,
                                        nest4_value,
                                    ) in nest4_data_dict.items():
                                        # find the last level attributes
                                        nest4_attributes[nest4_key] = attr.ib(
                                            default=nest4_value
                                        )

                                    #### LEVEL 4
                                    # 1. Create class,
                                    Nest4Class = attr.make_class(
                                        nest3_key, nest4_attributes, frozen=True
                                    )
                                    # 2. Update the data_dict to include the instantiated classes
                                    # 3. Instantiate the class
                                    nest4_instance = Nest4Class(**nest4_data_dict)
                                    # 4. Add the instantiated class to the dict
                                    nest4_instances[nest3_key] = nest4_instance
                                    # 5. Add to the attributes
                                    nest3_attributes[nest3_key] = attr.ib(
                                        default=nest4_instance
                                    )
                                else:
                                    nest3_attributes[nest3_key] = attr.ib(
                                        default=nest3_value
                                    )

                            ### LEVEL 3
                            # 1. Create class,
                            Nest3Class = attr.make_class(
                                nest2_key, nest3_attributes, frozen=True
                            )
                            # 2. Update the data_dict to include the instantiated classes
                            nest3_data_dict = update_dict_with_instantiated_classes(
                                nest3_data_dict, nest4_instances
                            )
                            # 3. Instantiate the class
                            nest3_instance = Nest3Class(**nest3_data_dict)
                            # 4. Add the instantiated class to the dict
                            nest3_instances[nest2_key] = nest3_instance
                            # 5. Add to the attributes
                            nest2_attributes[nest2_key] = attr.ib(
                                default=nest3_instance
                            )
                        else:
                            nest2_attributes[nest2_key] = attr.ib(default=nest2_value)

                    ### LEVEL 2
                    # 1. Create class,
                    Nest2Class = attr.make_class(
                        nest1_key, nest2_attributes, frozen=True
                    )
                    # 2. Update the data_dict to include the instantiated classes
                    nest2_data_dict = update_dict_with_instantiated_classes(
                        nest2_data_dict, nest3_instances
                    )
                    # 3. Instantiate the class
                    nest2_instance = Nest2Class(**nest2_data_dict)
                    # 4. Add the instantiated class to the dict
                    nest2_instances[nest1_key] = nest2_instance
                    # 5. Add to the attributes
                    nest1_attributes[nest1_key] = attr.ib(default=nest2_instance)
                else:
                    nest1_attributes[nest1_key] = attr.ib(default=nest1_value)

            ### LEVEL 1
            # 1. Create class,
            Nest1Class = attr.make_class(nest0_key, nest1_attributes, frozen=True)
            # 2. Update the data_dict to include the instantiated classes
            nest1_data_dict = update_dict_with_instantiated_classes(
                nest1_data_dict, nest2_instances
            )
            # 3. Instantiate the class
            nest1_instance = Nest1Class(**nest1_data_dict)
            # 4. Add the instantiated class to the dict
            nest1_instances[nest0_key] = nest1_instance
            # Add to the attributes
            nest0_attributes[nest0_key] = attr.ib(default=nest1_instance)
        else:
            nest0_attributes[nest0_key] = attr.ib(default=nest0_value)

    ### LEVEL 0
    # 1. Create class,
    Nest0Class = attr.make_class(class_name, nest0_attributes, frozen=True)
    # 2. Update the data_dict to include the instantiated classes
    nest0_data_dict = update_dict_with_instantiated_classes(
        nest0_data_dict, nest1_instances
    )
    # 3. Instantiate the class
    nest0_instance = Nest0Class(**nest0_data_dict)

    return nest0_instance
