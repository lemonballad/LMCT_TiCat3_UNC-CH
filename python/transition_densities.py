import re
import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial

# Paths and constants
READ_FLAG = False
WRITE_FLAG = False
PATH_READ = r'C:\Users\Thomas\Desktop\TiCat3\Gaussian\Good Result\El_C3'
PATH_WRITE = r'C:\Users\Thomas\Desktop\TiCat3\Gaussian\Good Result\El_C3'
MOLECULE = 'TiCat3_TD'
FILE_EXT_CUB = '.cub'
FILE_EXT_LOG = '.log'
CUBE_1 = f"{PATH_READ}\\{MOLECULE}_1{FILE_EXT_CUB}"
LOG_1 = f"{PATH_READ}\\{MOLECULE}_1{FILE_EXT_LOG}"


# Function to read cube data
def read_cube(current_cube):
    """
    Read data from a Gaussian cube file.

    Parameters:
    - current_cube (str): Path to the cube file.

    Returns:
    - _num_atoms_ (int): Number of atoms.
    - _origin_ (np.ndarray): Origin coordinates.
    - _dimensions_ (np.ndarray): Grid dimensions.
    - _scaling_ (np.ndarray): Grid scaling.
    - _atomic_numbers_ (np.ndarray): Atomic numbers.
    - _charges_ (np.ndarray): Atomic charges.
    - _positions_ (np.ndarray): Atomic positions.
    - density (np.ndarray): Electron density grid.
    """
    data = np.genfromtxt(current_cube, skip_header=1)
    _num_atoms_ = int(data[0, 0])
    _origin_ = data[0, 1:4]
    _dimensions_ = data[1:4, 0].astype(int)
    _scaling_ = data[1:4, 1:4]
    _atomic_numbers_ = data[5:5 + _num_atoms_, 0].astype(int)
    _charges_ = data[5:5 + _num_atoms_, 1]
    _positions_ = data[5:5 + _num_atoms_, 2:5]

    start_row = 6 + _num_atoms_
    temp = np.reshape(data[start_row:, 0], (_dimensions_[2], -1))
    next6 = _dimensions_[2] - (_dimensions_[2] % 6) + 6
    valid_indices = (np.arange(temp.size) % next6 < _dimensions_[2]) & (np.arange(temp.size) % next6 != 0)
    density = temp[valid_indices].reshape(_dimensions_[2], _dimensions_[1], _dimensions_[0])

    return _num_atoms_, _origin_, _dimensions_, _scaling_, _atomic_numbers_, _charges_, _positions_, density


# Function to write cube data
def write_cube(_header_, filename, density):
    """
    Write electron density to a Gaussian cube file.

    Parameters:
    - header (str): Header information.
    - filename (str): Output file path.
    - density (np.ndarray): Electron density grid.
    """
    with open(filename, 'w+') as file_out:
        file_out.write(_header_)
        for n1 in range(density.shape[2]):
            for n2 in range(density.shape[1]):
                for row in range(density.shape[0] // 6):
                    file_out.write(''.join([f'{density[row * 6 + i, n2, n1]:13.5E}' for i in range(6)]) + '\n')
                for row in range(density.shape[0] // 6 * 6, density.shape[0]):
                    file_out.write(
                        ''.join([f'{density[row, n2, n1]:13.5E}' for i in range(density.shape[0] % 6)]) + '\n')


# Regular expression pattern
EXCITATION_PATTERN = r' Excited State\s+(\d+):\s+\S+\s+(\d+\.\d+) eV\s+(\d+\.\d+) nm\s+f=(\d+\.\d+)'


# Function to parse log file
def parse_log_file(log_file):
    """
    Parse Gaussian log file to extract excitation energies, wavelengths, and oscillator strengths.

    Parameters:
    - log_file (str): Path to the log file.

    Returns:
    - _energies_ (np.ndarray): Excitation energies.
    - _wavelengths_ (np.ndarray): Corresponding wavelengths.
    - _oscillator_strengths_ (np.ndarray): Oscillator strengths.
    """
    _energies_ = np.zeros(20)
    _wavelengths_ = np.zeros(20)
    _oscillator_strengths_ = np.zeros(20)

    with open(log_file, 'r') as file_in:
        for line in file_in:
            match = re.search(EXCITATION_PATTERN, line)
            if match:
                _level_ = int(match.group(1))
                _energies_[_level_ - 1] = float(match.group(2))
                _wavelengths_[_level_ - 1] = float(match.group(3))
                _oscillator_strengths_[_level_ - 1] = float(match.group(4))
                if _level_ == 20:
                    break

    return _energies_, _wavelengths_, _oscillator_strengths_


# Function to condense transition density
def condense_transition_density(_center_, _radius_, _x_, _y_, _z_, rho):
    """
    Calculate condensed transition density and error term.

    Parameters:
    - _center_ (np.ndarray): Centers of the regions.
    - _radius_ (np.ndarray): Radii of the regions.
    - _x_ (np.ndarray): X-axis coordinates.
    - _y_ (np.ndarray): Y-axis coordinates.
    - _z_ (np.ndarray): Z-axis coordinates.
    - rho (np.ndarray): Transition density.

    Returns:
    - _cond_rho_ (np.ndarray): Condensed transition densities.
    - _err_term_ (float): Error term.
    """
    d_volume = abs(_x_[1] - _x_[0]) * abs(_y_[1] - _y_[0]) * abs(_z_[1] - _z_[0])
    _num_centers_ = len(_radius_)

    _cond_rho_ = np.zeros(4)
    remains = rho.copy()

    for ic in range(_num_centers_):
        y_mat, z_mat, x_mat = np.meshgrid((_y_ - _center_[1, ic]) ** 2, (_z_ - _center_[2, ic]) ** 2, (_x_ - _center_[0, ic]) ** 2)

        if ic == 0:
            dist = np.sqrt(x_mat + y_mat) - _radius_[ic]
        else:
            dist = np.sqrt(x_mat + y_mat + z_mat) - _radius_[ic]

        _cond_rho_[ic] = np.trapz(rho[dist < 0].flatten(), dx=d_volume)
        remains[dist < 0] = 0

    _err_term_ = np.trapz(np.trapz(np.trapz(remains, axis=0), axis=1), axis=2) * d_volume

    return _cond_rho_, _err_term_


# Parallel processing function for condensing transition density
def parallel_condense(args, _center_, _radius_, _x_, _y_, _z_):
    """
    Wrapper function for parallelizing condense_transition_density.

    Parameters:
    - _args_ (tuple): Tuple containing level number and transition density.
    - _center_ (np.ndarray): Centers of the regions.
    - _radius_ (np.ndarray): Radii of the regions.
    - _x_ (np.ndarray): X-axis coordinates.
    - _y_ (np.ndarray): Y-axis coordinates.
    - _z_ (np.ndarray): Z-axis coordinates.

    Returns:
    - _level_ (int): Energy level number.
    - _cond_rho_ (np.ndarray): Condensed transition densities.
    - _err_term_ (float): Error term.
    """
    _level_, rho = args
    return _level_, condense_transition_density(_center_, _radius_, _x_, _y_, _z_, rho)


# Main processing
if __name__ == '__main__':
    # Read cube data
    num_atoms, origin, dimensions, scaling, atomic_numbers, charges, positions, density_1 = read_cube(CUBE_1)

    # Read header lines from cube file
    with open(CUBE_1, 'r') as file_in:
        header = ''.join([next(file_in) for _ in range(6 + num_atoms)])

    # Parse log file
    energies, wavelengths, oscillator_strengths = parse_log_file(LOG_1)

    # Define center and radius arrays
    center = np.vstack((positions[0], np.mean(positions[[8, 11], :], axis=0), np.mean(positions[[19, 22], :], axis=0),
                        np.mean(positions[[29, 32], :], axis=0)))
    radius = np.array([4.0, 3.5, 3.5, 3.5])

    # Define grid coordinates
    X = scaling[0, 0] * np.arange(dimensions[0]) + origin[0]
    Y = scaling[1, 1] * np.arange(dimensions[1]) + origin[1]
    Z = scaling[2, 2] * np.arange(dimensions[2]) + origin[2]

    # Initialize arrays for transition densities
    if READ_FLAG:
        transition_densities = np.zeros((dimensions[2], dimensions[1], dimensions[0], 20), dtype=np.double)

    condensed_transition_densities = np.zeros((4, 20), dtype=np.double)
    error_terms = np.zeros(20)

    # Loop over electronic energy levels
    if cpu_count() > 1:  # Parallel processing if more than one CPU core available
        with Pool(processes=cpu_count()) as pool:
            results = pool.map(partial(parallel_condense, center=center, radius=radius, X=X, Y=Y, Z=Z),
                               enumerate(transition_densities.transpose(3, 2, 1, 0)))

        for result in results:
            level, (cond_rho, err_term) = result
            condensed_transition_densities[:, level] = cond_rho
            error_terms[level] = err_term
    else:
        for level in range(20):
            cond_rho, err_term = condense_transition_density(center, radius, X, Y, Z,
                                                             transition_densities[:, :, :, level])
            condensed_transition_densities[:, level] = cond_rho
            error_terms[level] = err_term

    # Write cube files if specified
    if WRITE_FLAG:
        for level in range(20):
            file_write = f"{PATH_WRITE}\\{MOLECULE}_{level + 1}_1{FILE_EXT_CUB}"
            write_cube(header, file_write, transition_densities[:, :, :, level])
