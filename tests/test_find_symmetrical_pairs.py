import numpy as np
from kitesim.structural.structural_mesher import find_symmetrical_pairs


## TODO: add tests with float working and without


def test_identifying_symmetrical_pairs():
    points = np.array(
        [
            [1, 1, 1],
            [1, -1, 1],
        ]
    )
    expected_value = np.array([[0, 1]])
    actual_value = find_symmetrical_pairs(points)
    assert isinstance(actual_value, np.ndarray)
    np.testing.assert_allclose(actual_value, expected_value)


def test_empty_array():
    points = np.array([])
    expected_value = np.array([])
    actual_value = find_symmetrical_pairs(points)
    assert isinstance(actual_value, np.ndarray)
    np.testing.assert_allclose(actual_value, expected_value)


def test_single_point():
    points = np.array([[1, 1, 1]])
    expected_value = np.array([])
    actual_value = find_symmetrical_pairs(points)
    assert isinstance(actual_value, np.ndarray)
    np.testing.assert_allclose(actual_value, expected_value)


def test_multiple_points():
    points = np.array(
        [
            [1, 1, 1],
            [2, 2, 2],
            [3, 3, 3],
            [4, 4, 4],
        ]
    )
    expected_value = np.array([])
    actual_value = find_symmetrical_pairs(points)
    assert isinstance(actual_value, np.ndarray)
    np.testing.assert_allclose(actual_value, expected_value)
