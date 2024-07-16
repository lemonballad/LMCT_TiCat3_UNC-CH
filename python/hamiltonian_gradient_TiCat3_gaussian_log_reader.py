import numpy as np
import re

# Constants
OMEGA_MAX = 25000
LN_WIDTH = 3000
FWHM = 2.355 * LN_WIDTH
DIST = 3.5848

# Regular expressions for parsing
E0_EXPR = r'SCF Done:\s+\S+\s+[=]+\s+([-]*\d+[.]+\d+)\s+'
MU_EXPR = r'\s+(?P<state>\d+)\s+([-]*\d+[.]+\d+)\s+([-]*\d+[.]+\d+)\s+([-]*\d+[.]+\d+)\s+\d+[.]+\d+\s+\d+[.]+\d+'
EX_STATE_EXPR = r'\s+Excited State\s+\d+[:]+\s+\S+\s+\d+[.]+\d+\s+eV\s+(?P<wln>\d+[.]+\d+)\s+nm\s+f[=]+(?P<Osc>\d+[.]+\d+)\s+'
EX_STATE_FLAG_EXPR = r' Excitation energies and oscillator strengths:'
EX_STATE_END_EXPR = r' Leave Link  914 '
MU_FLAG_EXPR = r'Ground to excited state transition electric dipole moments'
MU_END_EXPR = r'Ground to excited state transition velocity dipole moments'


def read_file(file_path):
    """Read a file and return its lines as a list.

    Args:
        file_path: Path to the file.

    Returns:
        list: Lines of the file.
    """
    with open(file_path, 'r') as file:
        return file.readlines()


def parse_energy(line):
    """Parse ground state energy from a line.

    Args:
        line: Line from the file.

    Returns:
        float: Ground state energy in a.u. or None if no match.
    """
    match = re.search(E0_EXPR, line)
    return 2.2 * 10 ** 5 * float(match.group(1)) if match else None


def parse_transition_dipole(line):
    """Parse transition dipole from a line.

    Args:
        line: Line from the file.

    Returns:
        tuple: State and dipole value, or (None, None) if no match.
    """
    match = re.search(MU_EXPR, line)
    if match:
        state = int(match.group('state'))
        mu_x, mu_y, mu_z = map(float, [match.group(2), match.group(3), match.group(4)])
        mu_val = np.sqrt(mu_x ** 2 + mu_y ** 2 + mu_z ** 2)
        return state, mu_val
    return None, None


def parse_excited_state(line):
    """Parse excited state energies and oscillator strengths from a line.

    Args:
        line: Line from the file.

    Returns:
        tuple: Wavelength and oscillator strength, or (None, None) if no match.
    """
    match = re.search(EX_STATE_EXPR, line)
    if match:
        wln = float(match.group('wln'))
        osc = float(match.group('Osc'))
        return wln, osc
    return None, None


def adjust_variable_length_data(data, default_val=np.nan):
    """Adjust for variable length data.

    Args:
        data (np.ndarray): Input data array.
        default_val: Default value for padding.

    Returns:
        np.ndarray: Adjusted data array.
    """
    data[1:, 1:4, :] = data[1:, 0:3, :]
    data[1:, [0, 4], :] = default_val
    return data


def compute_extinction_coefficient(dipole_moments):
    """Compute extinction coefficient.

    Args:
        dipole_moments: Transition dipole moments array.

    Returns:
        np.ndarray: Extinction coefficients array.
    """
    return dipole_moments ** 2 * OMEGA_MAX / 0.009186 / np.sqrt(2 * np.pi) / LN_WIDTH


def compute_electronic_coupling(extinction_coefficient):
    """Compute electronic coupling.

    Args:
        extinction_coefficient: Extinction coefficients array.

    Returns:
        np.ndarray: Electronic coupling array.
    """
    return 0.0206 / DIST * np.sqrt(extinction_coefficient * OMEGA_MAX * FWHM)


def compute_gradients(couplings, energies):
    """Compute gradients of couplings and energies.

    Args:
        couplings: Electronic coupling array.
        energies: Excited state energies array.

    Returns:
        tuple: Gradients of couplings and energies.
    """
    _dH_AD_ = np.abs(np.diff(couplings[:, 1:4, :], axis=1)) / 0.02
    _grad_E_ = np.abs(np.diff(energies[:, 1:4, :], axis=1)) / 0.02
    _dH_AD_mean_ = np.mean(_dH_AD_, axis=1)
    _grad_E_mean_ = np.mean(_grad_E_, axis=1)
    return _dH_AD_mean_, _grad_E_mean_


def main():
    """Main function to process Gaussian log files."""
    # Paths and file settings
    path_read = 'C:/Users/Thomas/Desktop/TiCat3/Gaussian/Output/3-17-16'
    molecule_prefix = 'TiCat3_1_'
    file_extension_log = '.log'

    # Normal modes
    normal_modes = [6, 10, 27, 39]

    # Initialize data matrices
    delta_excitation_energy = np.zeros((len(normal_modes), 5, 9), dtype=np.double)
    ground_state_energy = np.zeros((len(normal_modes), 5), dtype=np.double)
    excited_state_energy = np.zeros((len(normal_modes), 5, 9), dtype=np.double)
    transition_dipole = np.zeros((len(normal_modes), 5, 9), dtype=np.double)
    oscillator_strength = np.zeros((len(normal_modes), 5, 9), dtype=np.double)

    # Loop over normal modes
    for mode_index, mode in enumerate(normal_modes):
        # Define array of displacements in normal modes
        displacements = np.arange(-0.02, 0.03, 0.01) if mode == 6 else np.arange(-0.01, 0.02, 0.01)

        # Loop over displacement in normal modes
        for disp_index, displacement in enumerate(displacements):
            level = 0
            log_file = f"{path_read}/{molecule_prefix}{mode}_{displacement:.2f}{file_extension_log}"

            lines = read_file(log_file)
            excited_state_flag = False
            transition_dipole_flag = False

            for line in lines:
                # Parse ground state energy
                energy = parse_energy(line)
                if energy:
                    ground_state_energy[mode_index, disp_index] = energy

                # Check for transition dipole flag
                if re.search(MU_FLAG_EXPR, line):
                    transition_dipole_flag = True

                # Parse transition dipole moments
                if transition_dipole_flag:
                    state, dipole = parse_transition_dipole(line)
                    if state is not None:
                        transition_dipole[mode_index, disp_index, state] = dipole
                    if re.search(MU_END_EXPR, line) or state == 9:
                        transition_dipole_flag = False

                # Check for excited state flag
                if re.search(EX_STATE_FLAG_EXPR, line):
                    excited_state_flag = True

                # Parse excited state energies and oscillator strengths
                if excited_state_flag:
                    wln, osc = parse_excited_state(line)
                    if wln and osc:
                        level += 1
                        excited_state_energy[mode_index, disp_index, level - 1] = 10 ** 7 / wln
                        delta_excitation_energy[mode_index, disp_index, level - 1] = \
                            excited_state_energy[mode_index, disp_index, level - 1] + ground_state_energy[
                                mode_index, disp_index]
                        oscillator_strength[mode_index, disp_index, level - 1] = osc
                    if re.search(EX_STATE_END_EXPR, line) or level == 9:
                        excited_state_flag = False

    # Adjust for variable length data
    delta_excitation_energy = adjust_variable_length_data(delta_excitation_energy)
    ground_state_energy = adjust_variable_length_data(ground_state_energy)
    excited_state_energy = adjust_variable_length_data(excited_state_energy)
    transition_dipole = adjust_variable_length_data(transition_dipole)
    oscillator_strength = adjust_variable_length_data(oscillator_strength)

    # Compute extinction coefficient
    extinction_coefficient = compute_extinction_coefficient(transition_dipole)

    # Compute electronic coupling
    electronic_coupling = compute_electronic_coupling(extinction_coefficient)
    H_ad_rev = (excited_state_energy[:, :, 5] - excited_state_energy[:, :, 3]) / 3
    diag_el = (excited_state_energy[:, :, 5] + 2 * excited_state_energy[:, :, 3]) / 3

    # Compute gradients
    dH_AD, grad_E = compute_gradients(electronic_coupling, excited_state_energy)

    # More gradient calculations
    dhad = (electronic_coupling[:, 1, :] - electronic_coupling[:, 3, :]) / -0.02
    grade = (excited_state_energy[:, 1, :] - excited_state_energy[:, 3, :]) / -0.02
    dhadrev = (H_ad_rev[:, 1] - H_ad_rev[:, 3]) / -0.02


if __name__ == "__main__":
    main()
