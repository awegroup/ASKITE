import numpy as np


def calculate_elongation(
    points_new,
    wing_rest_lengths,
    bridle_rest_lengths,
    config,
):
    # Initialise
    wing_ci = config.kite.connectivity.wing_ci
    wing_cj = config.kite.connectivity.wing_cj
    bridle_ci = config.kite.connectivity.bridle_ci
    bridle_cj = config.kite.connectivity.bridle_cj
    tube_line_indices = config.kite.connectivity.tube_line_indices

    # dL_kite_stretch_lst, dL_kite_slack_lst = np.zeros(len(ci_kite)), np.zeros(len(ci_kite))
    (
        dL_tube_stretch_lst,
        dL_canopy_stretch_lst,
        dL_bridle_stretch_lst,
        dL_wing_stretch_lst,
    ) = ([], [], [], [])
    dL_tube_slack_lst, dL_canopy_slack_lst, dL_bridle_slack_lst, dL_wing_slack_lst = (
        [],
        [],
        [],
        [],
    )
    dL_bridle, dL_wing = [], []
    for i, (ci, cj) in enumerate(zip(wing_ci, wing_cj)):
        sep_vec = points_new[ci] - points_new[cj]
        sep = np.linalg.norm(sep_vec)
        rest_length = wing_rest_lengths[i]
        dL_perc = 100.0 * ((sep - rest_length) / rest_length)

        dL_wing.append(dL_perc)

        if dL_perc > 0:
            dL_wing_stretch_lst.append(dL_perc)

            if i in tube_line_indices:
                dL_tube_stretch_lst.append(dL_perc)
            else:
                dL_canopy_stretch_lst.append(dL_perc)
        else:
            dL_wing_slack_lst.append(dL_perc)
            if i in tube_line_indices:
                dL_tube_slack_lst.append(dL_perc)
            else:
                dL_canopy_slack_lst.append(dL_perc)

    for i, (ci, cj) in enumerate(zip(bridle_ci, bridle_cj)):
        sep_vec = points_new[ci] - points_new[cj]
        sep = np.linalg.norm(sep_vec)
        rest_length = bridle_rest_lengths[i]
        dL_perc = 100.0 * ((sep - rest_length) / rest_length)

        dL_bridle.append(dL_perc)

        if dL_perc > 0:
            dL_bridle_stretch_lst.append(dL_perc)
        else:
            dL_bridle_slack_lst.append(dL_perc)

    return (
        [
            dL_tube_stretch_lst,
            dL_canopy_stretch_lst,
            dL_bridle_stretch_lst,
            dL_wing_stretch_lst,
        ],
        [
            dL_tube_slack_lst,
            dL_canopy_slack_lst,
            dL_bridle_slack_lst,
            dL_wing_slack_lst,
        ],
        np.concatenate((dL_bridle, dL_wing)),
    )
