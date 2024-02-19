import numpy as np
import matplotlib.pyplot as plt


def airfoil_coeffs(alpha, coeffs):
    cl = np.interp(alpha * 180 / np.pi, coeffs[:, 0], coeffs[:, 1])
    cd = np.interp(alpha * 180 / np.pi, coeffs[:, 0], coeffs[:, 2])
    cm = np.interp(alpha * 180 / np.pi, coeffs[:, 0], coeffs[:, 3])

    return cl, cd, cm


def LEI_airf_coeff(t, k, alpha):
    """
    ----------
    t : TYPE
        DESCRIPTION.
    k : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.

    Returns
    -------
    Cl : TYPE
        DESCRIPTION.
    Cd : TYPE
        DESCRIPTION.
    Cm : TYPE
        DESCRIPTION.

    """
    C20 = -0.008011
    C21 = -0.000336
    C22 = 0.000992
    C23 = 0.013936
    C24 = -0.003838
    C25 = -0.000161
    C26 = 0.001243
    C27 = -0.009288
    C28 = -0.002124
    C29 = 0.012267
    C30 = -0.002398
    C31 = -0.000274
    C32 = 0
    C33 = 0
    C34 = 0
    C35 = -3.371000
    C36 = 0.858039
    C37 = 0.141600
    C38 = 7.201140
    C39 = -0.676007
    C40 = 0.806629
    C41 = 0.170454
    C42 = -0.390563
    C43 = 0.101966
    C44 = 0.546094
    C45 = 0.022247
    C46 = -0.071462
    C47 = -0.006527
    C48 = 0.002733
    C49 = 0.000686
    C50 = 0.123685
    C51 = 0.143755
    C52 = 0.495159
    C53 = -0.105362
    C54 = 0.033468
    C55 = -0.284793
    C56 = -0.026199
    C57 = -0.024060
    C58 = 0.000559
    C59 = -1.787703
    C60 = 0.352443
    C61 = -0.839323
    C62 = 0.137932

    S9 = C20 * t**2 + C21 * t + C22
    S10 = C23 * t**2 + C24 * t + C25
    S11 = C26 * t**2 + C27 * t + C28
    S12 = C29 * t**2 + C30 * t + C31
    S13 = C32 * t**2 + C33 * t + C34
    S14 = C35 * t**2 + C36 * t + C37
    S15 = C38 * t**2 + C39 * t + C40
    S16 = C41 * t**2 + C42 * t + C43

    lambda5 = S9 * k + S10
    lambda6 = S11 * k + S12
    lambda7 = S13 * k + S14
    lambda8 = S15 * k + S16

    Cl = lambda5 * alpha**3 + lambda6 * alpha**2 + lambda7 * alpha + lambda8
    Cd = (
        ((C44 * t + C45) * k**2 + (C46 * t + C47) * k + (C48 * t + C49)) * alpha**2
        + (C50 * t + C51) * k
        + (C52 * t**2 + C53 * t + C54)
    )
    Cm = (
        ((C55 * t + C56) * k + (C57 * t + C58)) * alpha**2
        + (C59 * t + C60) * k
        + (C61 * t + C62)
    )

    if alpha > 20 or alpha < -20:
        Cl = 2 * np.cos(alpha * np.pi / 180) * np.sin(alpha * np.pi / 180) ** 2
        Cd = 2 * np.sin(alpha * np.pi / 180) ** 3

    return Cl, Cd, Cm


def calculate_polar_lookup_table(
    controlpoints,
    n_segments,
    n_panels_aero,
    n_split,
    tube_diameters,
    is_tube_diameter_dimensionless,
    canopy_max_heights,
    is_canopy_max_height_dimensionless,
):
    """Returns the 2D coefficients of the kite
    using the interpolated correlations from [Breukels2011]

    input:
        controlpoints       : list of control points of each refined segment
        N_segments_refined  : number of refined segments (aero sections)
        N_segments          : number of segments (canopy pieces between struts
        n_split             : number of times each segment is splitted for the aero-refinement

    output:
        data_airf           : array [[alpha,Cl,Cd,Cm],..] of 2D coefficients for each refined segment
    """

    thicc = np.array([])
    for i in range(n_segments):
        temp = np.linspace(tube_diameters[i], tube_diameters[i + 1], n_split + 1)
        temp1 = []
        for a in range(len(temp) - 1):
            temp1 = np.append(temp1, (temp[a] + temp[a + 1]) / 2)
            # temp1.append((temp[a] +temp[a+1])/2)
        thicc = np.append(thicc, temp1)

    camb = np.array([])
    for i in range(n_segments):
        temp = np.linspace(
            canopy_max_heights[i], canopy_max_heights[i + 1], n_split + 1
        )
        temp1 = []
        for a in range(len(temp) - 1):
            temp1 = np.append(temp1, (temp[a] + temp[a + 1]) / 2)
            # temp1.append((temp[a] +temp[a+1])/2)
        camb = np.append(camb, temp1)

    aoas = np.arange(-20, 21, 1)
    data_airf = np.empty((len(aoas), 4, n_panels_aero))
    thicc_c = np.empty(n_panels_aero)
    camb_c = np.empty(n_panels_aero)
    for i in range(n_panels_aero):
        for j in range(len(aoas)):
            if is_tube_diameter_dimensionless:
                thicc_c[i] = thicc[i]
            else:
                thicc_c[i] = thicc[i] / controlpoints[i]["chord"]

            if is_canopy_max_height_dimensionless:
                camb_c[i] = camb[i]
            else:
                camb_c[i] = camb[i] / controlpoints[i]["chord"]

            alpha = aoas[j]
            Cl, Cd, Cm = LEI_airf_coeff(thicc_c[i], camb_c[i], alpha)
            data_airf[j, 0, i] = alpha
            data_airf[j, 1, i] = Cl
            data_airf[j, 2, i] = Cd
            data_airf[j, 3, i] = Cm

    return data_airf


def create_geometry_LEI(coordinates, Uinf, N, ring_geo, model):
    filaments = []
    controlpoints = []
    rings = []
    wingpanels = []
    ringvec = []
    coord_L = []
    for i in range(N - 1):
        section = {
            "p1": coordinates[2 * i, :],
            "p2": coordinates[2 * i + 2, :],
            "p3": coordinates[2 * i + 3, :],
            "p4": coordinates[2 * i + 1, :],
        }
        wingpanels.append(section)
        di = vec_norm(
            coordinates[2 * i, :] * 0.75
            + coordinates[2 * i + 1, :] * 0.25
            - (coordinates[2 * i + 2, :] * 0.75 + coordinates[2 * i + 3, :] * 0.25)
        )
        if i == 0:
            diplus = vec_norm(
                coordinates[2 * (i + 1), :] * 0.75
                + coordinates[2 * (i + 1) + 1, :] * 0.25
                - (
                    coordinates[2 * (i + 1) + 2, :] * 0.75
                    + coordinates[2 * (i + 1) + 3, :] * 0.25
                )
            )
            ncp = di / (di + diplus)
        elif i == N - 2:
            dimin = vec_norm(
                coordinates[2 * (i - 1), :] * 0.75
                + coordinates[2 * (i - 1) + 1, :] * 0.25
                - (
                    coordinates[2 * (i - 1) + 2, :] * 0.75
                    + coordinates[2 * (i - 1) + 3, :] * 0.25
                )
            )
            ncp = dimin / (dimin + di)
        else:
            dimin = vec_norm(
                coordinates[2 * (i - 1), :] * 0.75
                + coordinates[2 * (i - 1) + 1, :] * 0.25
                - (
                    coordinates[2 * (i - 1) + 2, :] * 0.75
                    + coordinates[2 * (i - 1) + 3, :] * 0.25
                )
            )
            diplus = vec_norm(
                coordinates[2 * (i + 1), :] * 0.75
                + coordinates[2 * (i + 1) + 1, :] * 0.25
                - (
                    coordinates[2 * (i + 1) + 2, :] * 0.75
                    + coordinates[2 * (i + 1) + 3, :] * 0.25
                )
            )
            ncp = 0.25 * (dimin / (dimin + di) + di / (di + diplus) + 1)

        ncp = 1 - ncp
        chord = np.linalg.norm(
            (section["p2"] + section["p1"]) / 2 - (section["p3"] + section["p4"]) / 2
        )
        LLpoint = (section["p2"] * (1 - ncp) + section["p1"] * ncp) * 3 / 4 + (
            section["p3"] * (1 - ncp) + section["p4"] * ncp
        ) * 1 / 4
        VSMpoint = (section["p2"] * (1 - ncp) + section["p1"] * ncp) * 1 / 4 + (
            section["p3"] * (1 - ncp) + section["p4"] * ncp
        ) * 3 / 4
        coord_L.append(LLpoint)

        # Define bound vortex filament
        bound = {
            "id": "bound",
            "x1": section["p1"] * 3 / 4 + section["p4"] * 1 / 4,
            "x2": section["p2"] * 3 / 4 + section["p3"] * 1 / 4,
            "Gamma": 0,
        }
        filaments.append(bound)

        x_airf = np.cross(VSMpoint - LLpoint, section["p2"] - section["p1"])
        x_airf = x_airf / np.linalg.norm(x_airf)
        y_airf = VSMpoint - LLpoint
        y_airf = y_airf / np.linalg.norm(y_airf)
        z_airf = bound["x2"] - bound["x1"]
        # z_airf[0] = 0
        z_airf = z_airf / np.linalg.norm(z_airf)
        airf_coord = np.column_stack([x_airf, y_airf, z_airf])

        normal = x_airf
        tangential = y_airf
        if model == "VSM":
            cp = {
                "coordinates": VSMpoint,
                "chord": chord,
                "normal": normal,
                "tangential": tangential,
                "airf_coord": airf_coord,
                "coordinates_aoa": LLpoint,
            }
            controlpoints.append(cp)
        elif model == "LLT":
            cp = {
                "coordinates": LLpoint,
                "chord": chord,
                "normal": normal,
                "tangential": tangential,
                "airf_coord": airf_coord,
            }
            controlpoints.append(cp)

        temp = {
            "r0": bound["x2"] - bound["x1"],
            "r1": cp["coordinates"] - bound["x1"],
            "r2": cp["coordinates"] - bound["x2"],
            "r3": cp["coordinates"] - (bound["x2"] + bound["x1"]) / 2,
        }
        ringvec.append(temp)

        temp = Uinf / np.linalg.norm(Uinf)
        if ring_geo == "3fil":
            # create trailing filaments, at x1 of bound filament
            temp1 = {"dir": temp, "id": "trailing_inf1", "x1": bound["x1"], "Gamma": 0}
            filaments.append(temp1)

            # create trailing filaments, at x2 of bound filament
            temp1 = {"x1": bound["x2"], "dir": temp, "id": "trailing_inf2", "Gamma": 0}
            filaments.append(temp1)
        elif ring_geo == "5fil":
            temp1 = {
                "x1": section["p4"],
                "x2": bound["x1"],
                "Gamma": 0,
                "id": "trailing1",
            }
            filaments.append(temp1)

            temp1 = {
                "dir": temp,
                "id": "trailing_inf1",
                "x1": section["p4"],
                "Gamma": 0,
            }
            filaments.append(temp1)

            # create trailing filaments, at x2 of bound filament
            temp1 = {
                "x2": section["p3"],
                "x1": bound["x2"],
                "Gamma": 0,
                "id": "trailing1",
            }
            filaments.append(temp1)

            temp1 = {
                "x1": section["p3"],
                "dir": temp,
                "id": "trailing_inf2",
                "Gamma": 0,
            }
            filaments.append(temp1)

        #

        rings.append(filaments)
        filaments = []

    coord_L = np.array(coord_L)
    return controlpoints, rings, wingpanels, ringvec, coord_L


def cross_product(r1, r2):
    """
    Cross product between r1 and r2

    """

    return np.array(
        [
            r1[1] * r2[2] - r1[2] * r2[1],
            r1[2] * r2[0] - r1[0] * r2[2],
            r1[0] * r2[1] - r1[1] * r2[0],
        ]
    )


def vec_norm(v):
    """
    Norm of a vector

    """
    return np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def refine_LEI_mesh_ballooning(coord, ball_angle, n_split):
    wingpanels = []
    for i in range(len(ball_angle)):
        section = {
            "p1": coord[2 * i, :],
            "p2": coord[2 * i + 2, :],
            "p3": coord[2 * i + 3, :],
            "p4": coord[2 * i + 1, :],
        }
        wingpanels.append(section)

    refined_coord = []
    for i_sec in range(len(wingpanels)):
        angle = ball_angle[i_sec] * np.pi / 180
        L_sec1 = vec_norm(wingpanels[i_sec]["p2"] - wingpanels[i_sec]["p1"])
        R1 = L_sec1 / 2 / np.sin(angle)
        L_sec2 = vec_norm(wingpanels[i_sec]["p3"] - wingpanels[i_sec]["p4"])
        R2 = L_sec2 / 2 / np.sin(angle)
        zvec = (wingpanels[i_sec]["p2"] + wingpanels[i_sec]["p1"]) / 2 - (
            wingpanels[i_sec]["p4"] + wingpanels[i_sec]["p3"]
        ) / 2
        zvec = zvec / vec_norm(zvec)

        xvec1 = wingpanels[i_sec]["p2"] - wingpanels[i_sec]["p1"]
        xvec1 = xvec1 / vec_norm(xvec1)
        yvec1 = cross_product(zvec, xvec1)
        yvec1 = yvec1 / vec_norm(yvec1)

        xvec2 = wingpanels[i_sec]["p3"] - wingpanels[i_sec]["p4"]
        xvec2 = xvec2 / vec_norm(xvec2)
        yvec2 = cross_product(zvec, xvec2)
        yvec2 = yvec2 / vec_norm(yvec2)

        if i_sec > 4:
            xvec1 = wingpanels[i_sec]["p1"] - wingpanels[i_sec]["p2"]
            xvec1 = xvec1 / vec_norm(xvec1)
            xvec2 = wingpanels[i_sec]["p4"] - wingpanels[i_sec]["p3"]
            xvec2 = xvec2 / vec_norm(xvec2)

        xloc1 = np.linspace(-L_sec1 / 2, L_sec1 / 2, n_split)
        yloc01 = np.sqrt(R1**2 - (L_sec1 / 2) ** 2)
        yloc1 = -np.sqrt(R1**2 - xloc1**2) + yloc01
        zloc1 = np.zeros(n_split)

        xloc2 = np.linspace(-L_sec2 / 2, L_sec2 / 2, n_split)
        yloc02 = np.sqrt(R2**2 - (L_sec2 / 2) ** 2)
        yloc2 = -np.sqrt(R2**2 - xloc2**2) + yloc02
        zloc2 = np.zeros(n_split)

        vec1 = np.array([xvec1, yvec1, zvec]).T
        vec2 = np.array([xvec2, yvec2, zvec]).T
        ax_pos1 = (wingpanels[i_sec]["p2"] + wingpanels[i_sec]["p1"]) / 2
        ax_pos2 = (wingpanels[i_sec]["p3"] + wingpanels[i_sec]["p4"]) / 2
        temp_coord = np.empty((int(n_split * 2), 3))
        for i_spl in range(n_split):
            coord_loc1 = np.array([xloc1[i_spl], yloc1[i_spl], zloc1[i_spl]])
            coord_loc2 = np.array([xloc2[i_spl], yloc2[i_spl], zloc2[i_spl]])
            coord1 = np.matmul(vec1, coord_loc1) + ax_pos1
            coord2 = np.matmul(vec2, coord_loc2) + ax_pos2

            if i_sec > 4:
                ind = 2 * n_split - 1 - (2 * i_spl + 1)
                temp_coord[ind] = coord1
                ind = 2 * n_split - 1 - 2 * i_spl
                temp_coord[ind] = coord2
            else:
                temp_coord[2 * i_spl] = coord1
                temp_coord[2 * i_spl + 1] = coord2

        if i_sec == 0:
            refined_coord = temp_coord
        else:
            refined_coord = np.append(refined_coord, temp_coord[2::, :], axis=0)

    return refined_coord


def plot_polars(data_airf, section_number):
    # Extract alpha values
    alpha_values = data_airf[:, 0, section_number]

    # Extract cl, cd, cm values
    cl_values = data_airf[:, 1, section_number]
    cd_values = data_airf[:, 2, section_number]
    cm_values = data_airf[:, 3, section_number]

    # print(f"Alpha values: {alpha_values}")
    # print(f"Cl values: {cl_values}")
    # print(f"Cd values: {cd_values}")
    # print(f"Cm values: {cm_values}")

    # Plot cl-alpha
    plt.figure(figsize=(10, 6))
    for i in range(len(alpha_values)):
        plt.scatter(
            alpha_values[i],
            cl_values[i],
            label=f"Alpha = {alpha_values[i]}",
            color="red",
        )
    plt.xlabel("Alpha")
    plt.ylabel("Cl")
    plt.title("Cl vs Alpha")
    # plt.legend()
    plt.grid(True)
    plt.show()

    # Plot cd-alpha
    plt.figure(figsize=(10, 6))
    for i in range(len(alpha_values)):
        plt.scatter(
            alpha_values[i],
            cd_values[i],
            label=f"Alpha = {alpha_values[i]}",
            color="red",
        )
    plt.xlabel("Alpha")
    plt.ylabel("Cd")
    plt.title("Cd vs Alpha")
    # plt.legend()
    plt.grid(True)
    plt.show()

    # Plot cm-alpha
    plt.figure(figsize=(10, 6))
    for i in range(len(alpha_values)):
        plt.scatter(
            alpha_values[i],
            cm_values[i],
            label=f"Alpha = {alpha_values[i]}",
            color="red",
        )
    plt.xlabel("Alpha")
    plt.ylabel("Cm")
    plt.title("Cm vs Alpha")
    # plt.legend()
    plt.grid(True)
    plt.show()
