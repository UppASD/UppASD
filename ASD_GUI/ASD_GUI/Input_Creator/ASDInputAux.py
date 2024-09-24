"""@package ASDInputAux
Set of auxiliary functions to write restartfiles .
It has functions to generate the coordinate file, as well as to write the
following magnetic configurations:
    - Domain walls:
        - Neel planar wall.
        - Bloch planar wall.
        - Vortex wall.
    - Skyrmion states:
        - Neel skyrmion.
        - Bloch skyrmion.
    - Helical spin spirals.

Author
----------
Jonathan Chico
"""
import numpy as np

##########################################################################
# @brief Function to generate domain wall profiles.
#
# @details It can generate
# - Planar domain walls:
#     - Neel planar wall.
#     - Bloch planar wall.
# - Vortex domain walls.
#
# The domain walls can be centered at a certain point in the lattice and with a given width.
# The chirality and type of wall (Neel and Bloch) can be chosen.
# @author Jonathan Chico
##########################################################################


def write_domain_wall(Natom, Mensemble, coord, DWInfo):
    """
    Generate magnetization for a domain wall in a system of atoms.

    Parameters:
    Natom (int): Number of atoms.
    Mensemble (int): Number of ensembles.
    coord (ndarray): Atomic coordinates.
    DWInfo (dict): Domain wall information.

    Returns:
    ndarray: Magnetization array.
    """

    tol = 1e-9
    # Magnetization of the system
    mag = np.zeros([Mensemble, Natom, 3], dtype=np.float64)
    # Setup a planar domain wall type, i.e. a 1D object
    if DWInfo["type"] == "planar":
        # Plane indicates which is the "propagation direction of the DW"
        # The sign in front of the 1/cosh determines the chirality of the wall
        if DWInfo["DWtype"] == "Bloch":
            if DWInfo["easy_axis"] == "x":
                # ---------------------------------------------------------------
                # Loop over the atoms in the system
                # ---------------------------------------------------------------
                for jj in range(0, Mensemble):
                    for ii in range(0, Natom):
                        arg = (DWInfo["center"] - coord[ii, DWInfo["plane"]]) / DWInfo[
                            "width"
                        ]
                        mag[:, ii, 0] = np.tanh(arg)
                        mag[:, ii, 1] = 0.0
                        mag[:, ii, 2] = DWInfo["chirality"] * 1.0 / np.cosh(arg)
                        mod = np.sqrt(mag[jj, ii].dot(mag[jj, ii]))
                        # Normalization of the spins
                        mag[jj, ii, :] = mag[jj, ii, :] / mod
            elif DWInfo["easy_axis"] == "y":
                # ---------------------------------------------------------------
                # Loop over the atoms in the system
                # ---------------------------------------------------------------
                for jj in range(0, Mensemble):
                    for ii in range(0, Natom):
                        arg = (DWInfo["center"] - coord[ii, DWInfo["plane"]]) / DWInfo[
                            "width"
                        ]
                        mag[:, ii, 0] = 0.0
                        mag[:, ii, 1] = np.tanh(arg)
                        mag[:, ii, 2] = DWInfo["chirality"] * 1.0 / np.cosh(arg)
                        mod = np.sqrt(mag[jj, ii].dot(mag[jj, ii]))
                        # Normalization of the spins
                        mag[jj, ii, :] = mag[jj, ii, :] / mod
            elif DWInfo["easy_axis"] == "z":
                # ---------------------------------------------------------------
                # Loop over the atoms in the system
                # ---------------------------------------------------------------
                for jj in range(0, Mensemble):
                    for ii in range(0, Natom):
                        arg = (DWInfo["center"] - coord[ii, DWInfo["plane"]]) / DWInfo[
                            "width"
                        ]
                        mag[:, ii, 0] = 0.0
                        mag[:, ii, 1] = DWInfo["chirality"] * 1.0 / np.cosh(arg)
                        mag[:, ii, 2] = np.tanh(arg)
                        mod = np.sqrt(mag[jj, ii].dot(mag[jj, ii]))
                        # Normalization of the spins
                        mag[jj, ii, :] = mag[jj, ii, :] / mod
        elif DWInfo["DWtype"] == "Neel":
            if DWInfo["easy_axis"] == "x":
                # ---------------------------------------------------------------
                # Loop over the atoms in the system
                # ---------------------------------------------------------------
                for jj in range(0, Mensemble):
                    for ii in range(0, Natom):
                        arg = (DWInfo["center"] - coord[ii, DWInfo["plane"]]) / DWInfo[
                            "width"
                        ]
                        mag[:, ii, 0] = np.tanh(arg)
                        mag[:, ii, 1] = DWInfo["chirality"] * 1.0 / np.cosh(arg)
                        mag[:, ii, 2] = 0.0
                        mod = np.sqrt(mag[jj, ii].dot(mag[jj, ii]))
                        # Normalization of the spins
                        mag[jj, ii, :] = mag[jj, ii, :] / mod
            elif DWInfo["easy_axis"] == "y":
                # ---------------------------------------------------------------
                # Loop over the atoms in the system
                # ---------------------------------------------------------------
                for jj in range(0, Mensemble):
                    for ii in range(0, Natom):
                        arg = (DWInfo["center"] - coord[ii, DWInfo["plane"]]) / DWInfo[
                            "width"
                        ]
                        mag[jj, ii, 0] = DWInfo["chirality"] * 1.0 / np.cosh(arg)
                        mag[jj, ii, 1] = np.tanh(arg)
                        mag[jj, ii, 2] = 0
                        mod = np.sqrt(mag[jj, ii].dot(mag[jj, ii]))
                        # Normalization of the spins
                        mag[jj, ii, :] = mag[jj, ii, :] / mod
            elif DWInfo["easy_axis"] == "z":
                # ---------------------------------------------------------------
                # Loop over the atoms in the system
                # ---------------------------------------------------------------
                for jj in range(0, Mensemble):
                    for ii in range(0, Natom):
                        arg = (DWInfo["center"] - coord[ii, DWInfo["plane"]]) / DWInfo[
                            "width"
                        ]
                        mag[jj, ii, 0] = DWInfo["chirality"] * 1.0 / np.cosh(arg)
                        mag[jj, ii, 1] = 0.0
                        mag[jj, ii, 2] = np.tanh(arg)
                        mod = np.sqrt(mag[jj, ii].dot(mag[jj, ii]))
                        # Normalization of the spins
                        mag[jj, ii, :] = mag[jj, ii, :] / mod
    # Setup a vortex domain wall type, i.e. a 2D object
    elif DWInfo["type"] == "vortex":
        # For the vortex wall one must introduce a rotation like term
        # -----------------------------------------------------------------------
        # Loop over the atoms in the system
        # -----------------------------------------------------------------------
        for ii in range(0, Natom):
            r_x = coord[ii, 0] - DWInfo["center"][0]
            r_y = coord[ii, 1] - DWInfo["center"][1]
            mod_r2 = np.sqrt(r_x**2 + r_y**2)
            arg = mod_r2 / DWInfo["radius"]
            if mod_r2 > tol:
                theta = np.arctan2(r_x, r_y)
            else:
                theta = 0.0
            mag[:, ii, 0] = DWInfo["chirality"] * np.cos(theta) * np.tanh(arg)
            mag[:, ii, 1] = DWInfo["chirality"] * np.sin(theta) * np.tanh(arg)
            mag[:, ii, 2] = DWInfo["polarity"] * 1.0 / np.cosh(arg)
            mod = np.sqrt(mag[jj, ii].dot(mag[jj, ii]))
            # Normalization of the spins
            mag[jj, ii, :] = mag[jj, ii, :] / mod
    return mag


##########################################################################
# @brief Writes a restartfile with a skyrmion profile.
# @details Generates a skyrmion profile making use of the skyrmion profiles defined in
# Nat. Commun. 7, 13613 (2016):
#
# @f$m_x = \cos(m\phi+\gamma)\sin(\theta(r))@f$
#
# @f$m_y = \sin(m\phi+\gamma)\sin(\theta(r))@f$
#
# @f$m_z = \cos(\theta(r))@f$
#
# With the out of plane angle \f$\theta\f$ being given by
# @f$\theta(r)=\pi + \arcsin(\tanh((r+c)/w)) + \arcsin(\tanh((r-c)/w))@f$
#
# with @f$c@f$ being the center of the skyrmion and @f$w@f$ the radius.
# The polar angle @f$\phi@f$ is determined from the skyrmion center, @f$\gamma@f$
# determined the type of skyrmion with @f$\gamma=0@f$ implies a Neel Skyrmion and
# @f$\gamma=\frac{\pi}{2}@f$ is a Bloch Skyrmion.
# The skyrmion is assumed to have its core parallel to the z-axis.
#
# @author Jonathan Chico
##########################################################################


def write_skyrmion(Natom, Mensemble, coord, SkxInfo):
    """Generates a skyrmion profile making use of the skyrmion profiles defined in
    Nat. Commun. 7, 13613 (2016):
    * .math: m_x = \\cos(m\\phi+\\gamma)\\sin(\theta(r))
    * .math: m_y = \\sin(m\\phi+\\gamma)\\sin(\theta(r))
    * .math: m_z = \\cos(\theta(r))
    """

    tol = 1e-9
    # Magnetization of the system
    mag = np.zeros([Mensemble, Natom, 3], dtype=np.float64)
    # ---------------------------------------------------------------------------
    # Loop over the atoms in the system
    # ---------------------------------------------------------------------------
    for jj in range(0, Mensemble):
        for ii in range(0, Natom):
            # Define the factors for the in-plane phi angle
            r_x = coord[ii, 0] - SkxInfo["center"][0]
            r_y = coord[ii, 1] - SkxInfo["center"][1]
            mod_r = np.sqrt(r_x * r_x + r_y * r_y)
            # The distance from the center of the skyrmion to which the DW
            # walls are set
            rad = SkxInfo["width"] * 0.5
            r_p = (mod_r + rad) / SkxInfo["width"]
            r_n = (mod_r - rad) / SkxInfo["width"]
            # Out of plane angle
            theta = (
                np.pi
                + np.arcsin(np.tanh(r_p))
                + np.arcsin(np.tanh(r_n))
                + SkxInfo["polarity"]
            )
            # In-plane skyrmion profile angle
            if mod_r > tol:
                phi = np.arctan2(r_y, r_x)
            else:
                phi = 0.0
            # Magnetization per site
            mag[jj, ii, 0] = (
                SkxInfo["handness"]
                * np.sin(theta)
                * np.cos(SkxInfo["order"] * phi + SkxInfo["type"])
            )
            mag[jj, ii, 1] = (
                SkxInfo["handness"]
                * np.sin(theta)
                * np.sin(SkxInfo["order"] * phi + SkxInfo["type"])
            )
            mag[jj, ii, 2] = np.cos(theta)
            mod = np.sqrt(mag[jj, ii].dot(mag[jj, ii]))
            # Normalization of the spins
            mag[jj, ii, :] = mag[jj, ii, :] / mod
    return mag


##########################################################################
# @brief Function to generate a generalized spin spiral configuration
# @details This function creates a generalized spin spiral via Rodrigues rotations
# @f$ \mathbf{v}_{rot}= \mathbf{v}\cos\theta
# +\left(\mathbf{v}\times\mathbf{k}\right)\sin\theta
# +\mathbf{k}\left(\mathbf{k}\cdot\mathbf{v}\right)\left(1-\cos\theta\right)@f$
# First a rotation is done to generate the cone angle, this is done by generating an
# axis that is perpendicular to the \textbf{pitch vector} (that is @f$v@f$ in the Rodrigues formula)
# which is dubbed the \textbf{cone axis}.
# The pitch vector is then rotated by that axis by an angle @f$\theta@f$, i.e. the cone angle.
# The obtained vector is then the initial spin direction, that will be rotated by
# the pitch vector, with an angle given by @f$\theta= \mathbf{q}\cdot\mathbf{r}@f$
# with @f$\mathbf{q}@f$ being the spiral wavevector.
# @author Jonathan Chico
##########################################################################
# pylint: disable=invalid-name, no-name-in-module, no-member


def create_spiral(Natom, Mensemble, coord, HLInfo):
    """
    Create a helical spiral magnetic configuration.

    Args:
        Natom (int): Number of atoms.
        Mensemble (int): Number of ensembles.
        coord (np.ndarray): Atomic coordinates.
        HLInfo (dict): Helical spiral parameters including 'pitch_vector',
                       'prop_vector', 'cone_angle', and 'handness'.

    Returns:
        np.ndarray: Magnetic configuration array of shape (Mensemble, Natom, 3).
    """

    # Transform lists to np arrays to avoid problems
    rot_vector = np.asarray(HLInfo["pitch_vector"], dtype=np.float64)
    prop_vector = np.asarray(HLInfo["prop_vector"], dtype=np.float64)
    # Magnetization of the system
    mag = np.zeros([Mensemble, Natom, 3], dtype=np.float64)
    # ---------------------------------------------------------------------------
    # First do a rotation to find the rotates spin for the conical phase
    # ---------------------------------------------------------------------------
    # First create a vector perpendicular to the rotation axis
    test_r = np.random.rand(3)
    # Normalize the vector
    test_r = test_r / np.sqrt(test_r.dot(test_r))
    # Axis which one will use to rotate the spins to get the conical phase
    cone_axis = np.cross(rot_vector, test_r)
    # ---------------------------------------------------------------------------
    # Rotate the spin first to find the needed cone angle using Rodriges rotation
    # ---------------------------------------------------------------------------
    init_spin = (
        rot_vector * np.cos(HLInfo["cone_angle"])
        + np.cross(cone_axis, rot_vector) * np.sin(HLInfo["cone_angle"])
        + cone_axis * (cone_axis.dot(rot_vector)) * (1 - np.cos(HLInfo["cone_angle"]))
    )
    # ---------------------------------------------------------------------------
    # Loop over the ensembles and atoms of the system
    # ---------------------------------------------------------------------------
    for jj in range(0, Mensemble):
        for ii in range(0, Natom):
            #
            theta = prop_vector.dot(coord[ii, :]) * 2.0 * np.pi * HLInfo["handness"]
            # -------------------------------------------------------------------
            # Do a Rodrigues rotation to get the helical spiral state
            # -------------------------------------------------------------------
            mag[jj, ii, :] = (
                init_spin * np.cos(theta)
                + np.cross(rot_vector, init_spin) * np.sin(theta)
                + rot_vector * (rot_vector.dot(init_spin)) * (1 - np.cos(theta))
            )
            mod = np.sqrt(mag[jj, ii].dot(mag[jj, ii]))
            # Normalization of the spins
            mag[jj, ii, :] = mag[jj, ii, :] / mod
    return mag


##########################################################################
# @brief Function to generate the coordinated for a given lattice following the
# same structure than in \c UppASD.
# @details Routine taken from the \c UppASD \c geometry.f90
#
# @author Anders Bergman and Johan Hellsvik.
# @note Adapted to python by Jonathan Chico
##########################################################################


def create_coord(cell, ncell, Bas, block_size, mom):
    """
    Generate coordinates and magnetic moments for a crystal structure.

    Args:
        cell (np.ndarray): 3x3 matrix defining the unit cell vectors.
        ncell (tuple): Number of cells in each direction (nx, ny, nz).
        Bas (np.ndarray): Basis vectors of the atoms in the unit cell.
        block_size (int): Size of the blocks for iteration.
        mom (np.ndarray): Magnetic moments of the atoms in the unit cell.

    Returns:
        tuple: Coordinates and magnetic moments of the atoms in the crystal.
    """

    tol = 1e-9

    NA = len(Bas)

    detmatrix = (
        cell[0, 0] * cell[1, 1] * cell[2, 2]
        - cell[0, 0] * cell[1, 2] * cell[2, 1]
        + cell[0, 1] * cell[1, 2] * cell[2, 0]
        - cell[0, 1] * cell[1, 0] * cell[2, 2]
        + cell[0, 2] * cell[1, 0] * cell[2, 1]
        - cell[0, 2] * cell[1, 1] * cell[2, 0]
    )

    invmatrix = np.zeros([3, 3], dtype=np.float64)
    if abs(detmatrix) > tol:
        invmatrix[0, 0] = (
            cell[1, 1] * cell[2, 2] - cell[2, 1] * cell[1, 2]
        ) / detmatrix
        invmatrix[0, 1] = (
            cell[0, 2] * cell[2, 1] - cell[2, 2] * cell[0, 1]
        ) / detmatrix
        invmatrix[0, 2] = (
            cell[0, 1] * cell[1, 2] - cell[1, 1] * cell[0, 2]
        ) / detmatrix
        invmatrix[1, 0] = (
            cell[1, 2] * cell[2, 0] - cell[2, 2] * cell[1, 0]
        ) / detmatrix
        invmatrix[1, 1] = (
            cell[0, 0] * cell[2, 2] - cell[2, 0] * cell[0, 2]
        ) / detmatrix
        invmatrix[1, 2] = (
            cell[0, 2] * cell[1, 0] - cell[1, 2] * cell[0, 0]
        ) / detmatrix
        invmatrix[2, 0] = (
            cell[1, 0] * cell[2, 1] - cell[2, 0] * cell[1, 1]
        ) / detmatrix
        invmatrix[2, 1] = (
            cell[0, 1] * cell[2, 0] - cell[2, 1] * cell[0, 0]
        ) / detmatrix
        invmatrix[2, 2] = (
            cell[0, 0] * cell[1, 1] - cell[1, 0] * cell[0, 1]
        ) / detmatrix

    icvec = np.zeros([3], dtype=np.float64)
    bsf = np.zeros([3], dtype=np.float64)
    for I0 in range(0, NA):
        icvec[0] = (
            Bas[I0, 0] * invmatrix[0, 0]
            + Bas[I0, 1] * invmatrix[1, 0]
            + Bas[I0, 2] * invmatrix[2, 0]
        )
        icvec[1] = (
            Bas[I0, 0] * invmatrix[0, 1]
            + Bas[I0, 1] * invmatrix[1, 1]
            + Bas[I0, 2] * invmatrix[2, 1]
        )
        icvec[2] = (
            Bas[I0, 0] * invmatrix[0, 2]
            + Bas[I0, 1] * invmatrix[1, 2]
            + Bas[I0, 2] * invmatrix[2, 2]
        )
        bsf[0] = np.floor(icvec[0] + 1e-7)
        bsf[1] = np.floor(icvec[1] + 1e-7)
        bsf[2] = np.floor(icvec[2] + 1e-7)
        for mu in range(0, 3):
            Bas[I0, mu] = (
                Bas[I0, mu]
                - bsf[0] * cell[0, mu]
                - bsf[1] * cell[1, mu]
                - bsf[2] * cell[2, mu]
            )

    ii = 0
    coord = np.zeros([NA * ncell[0] * ncell[1] * ncell[2], 3], dtype=np.float64)
    mom_mag = np.zeros([NA * ncell[0] * ncell[1] * ncell[2]], dtype=np.float64)
    for II3 in range(0, ncell[2], block_size):
        for II2 in range(0, ncell[1], block_size):
            for II1 in range(0, ncell[0], block_size):
                for I3 in range(II3, min(II3 + block_size, ncell[2])):
                    for I2 in range(II2, min(II2 + block_size, ncell[1])):
                        for I1 in range(II1, min(II1 + block_size, ncell[0])):
                            for I0 in range(0, NA):
                                mom_mag[ii] = mom[I0]
                                for mu in range(0, 3):
                                    coord[ii, mu] = (
                                        I1 * cell[0, mu]
                                        + I2 * cell[1, mu]
                                        + I3 * cell[2, mu]
                                        + Bas[I0, mu]
                                    )
                                ii = ii + 1
    return coord, mom_mag
