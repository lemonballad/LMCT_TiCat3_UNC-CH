import os
import re
import shutil
import json


def read_config(config_path):
    """
    Reads the configuration file.

    Args:
        config_path (str): Path to the configuration file.

    Returns:
        dict: Configuration dictionary.
    """
    with open(config_path, 'r') as config_file:
        return json.load(config_file)


def read_normal_mode_displacements(log_file_path):
    """
    Reads the normal mode displacements from a Gaussian log file.

    Args:
        log_file_path (str): Path to the Gaussian log file.

    Returns:
        dict: Dictionary of normal mode coordinates.
        int: Number of atoms.
    """
    start_pattern = ' Full mass-weighted force constant matrix:'
    end_pattern = ' - Thermochemistry -'
    atom_line_pattern = r'\s+(?P<atom_index>\d+)\s+\d+\s+'
    displacement_patterns = [r'(?P<x>-*\d+\.\d+)\s+(?P<y>-*\d+\.\d+)\s+(?P<z>-*\d+\.\d+)'] * 3
    combined_pattern = atom_line_pattern + r'\s+'.join(displacement_patterns)

    mode_coordinates = {}
    reading_displacements = False
    last_atom_index = -1
    mode_index = 0

    try:
        with open(log_file_path, 'r') as log_file:
            for line in log_file:
                if re.search(start_pattern, line):
                    reading_displacements = True
                if reading_displacements and not re.search(end_pattern, line):
                    match_result = re.search(combined_pattern, line)
                    if match_result:
                        atom_index = int(match_result.group('atom_index'))
                        if atom_index not in mode_coordinates:
                            mode_coordinates[atom_index] = {}
                        for i in range(3):
                            displacement_index = mode_index + i + 1
                            mode_coordinates[atom_index][displacement_index] = [
                                float(match_result.group(f'x{i + 1}')),
                                float(match_result.group(f'y{i + 1}')),
                                float(match_result.group(f'z{i + 1}'))
                            ]
                        if atom_index < last_atom_index:
                            mode_index += 3
                        last_atom_index = atom_index
                if re.search(end_pattern, line):
                    break
    except FileNotFoundError:
        print(f"Error: The file {log_file_path} was not found.")
        return None, None

    return mode_coordinates, len(mode_coordinates)


def read_ground_state_coordinates(gjf_file_path, num_atoms):
    """
    Reads the ground state Cartesian coordinates from a Gaussian input file.

    Args:
        gjf_file_path (str): Path to the Gaussian input file.
        num_atoms (int): Number of atoms in the molecule.

    Returns:
        list: Atom labels.
        list: x coordinates.
        list: y coordinates.
        list: z coordinates.
        str: Header lines of the file.
        str: Tail lines of the file.
    """
    _atom_labels_ = [''] * num_atoms
    _x_coordinates_ = [0.0] * num_atoms
    _y_coordinates_ = [0.0] * num_atoms
    _z_coordinates_ = [0.0] * num_atoms
    atom_index = 0
    coordinate_pattern = r'\s+(?P<atom_label>\S+)\s+(?P<x>-*\d+\.\d+)\s+(?P<y>-*\d+\.\d+)\s+(?P<z>-*\d+\.\d+)'
    _header_lines_ = ''
    _tail_lines_ = ''
    reading_header = True

    try:
        with open(gjf_file_path, 'r') as gjf_file:
            for line in gjf_file:
                if len(line.strip()) == 58:
                    atom_index += 1
                    match_result = re.search(coordinate_pattern, line)
                    if match_result:
                        _atom_labels_[atom_index - 1] = match_result.group('atom_label')
                        _x_coordinates_[atom_index - 1] = float(match_result.group('x'))
                        _y_coordinates_[atom_index - 1] = float(match_result.group('y'))
                        _z_coordinates_[atom_index - 1] = float(match_result.group('z'))
                    reading_header = False
                elif reading_header:
                    _header_lines_ += line + '\n'
                else:
                    _tail_lines_ += line + '\n'
    except FileNotFoundError:
        print(f"Error: The file {gjf_file_path} was not found.")
        return None, None, None, None, None, None

    return (_atom_labels_, _x_coordinates_, _y_coordinates_, _z_coordinates_, _header_lines_,
            _tail_lines_)


def write_gjf_files(write_path, molecule, modes, _displacements_, atoms, x, y, z,
                    _mode_coordinates_, header, tail, chk_file):
    """
    Writes the Gaussian input files for different modes and displacements.

    Args:
        write_path (str): Path to write the Gaussian input files.
        molecule (str): Molecule name.
        modes (list): List of normal modes of interest.
        _displacements_ (list): List of displacements to apply.
        atoms (list): Atom labels.
        x (list): x coordinates.
        y (list): y coordinates.
        z (list): z coordinates.
        _mode_coordinates_ (dict): Normal mode coordinates.
        header (str): Header lines of the Gaussian input file.
        tail (str): Tail lines of the Gaussian input file.
        chk_file (str): Path to the checkpoint file.
    """
    for root in range(1, 2):  # 1 to 6 in original code, changed to 1 for simplicity
        for mode in modes:
            for displacement in _displacements_:
                header_current = re.sub(r'%chk=TiCat3_C3_0.chk',
                                        f'%chk=TiCat3_C3_{mode}_{displacement}.chk', header)

                body_current = ''
                for index in range(len(atoms)):
                    line = f' {atoms[index]:<16}{x[index] + displacement * _mode_coordinates_[index + 1][mode + 1][0]:14.8f}{y[index] + displacement * _mode_coordinates_[index + 1][mode + 1][1]:14.8f}{z[index] + displacement * _mode_coordinates_[index + 1][mode + 1][2]:14.8f}\n'
                    body_current += line

                file_content = header_current + body_current + tail

                gjf_file_path = os.path.join(write_path, f'{molecule}_{mode}_{displacement}.gjf')
                new_chk_file_path = os.path.join(write_path, f'{molecule}_{mode}_{displacement}.chk')

                try:
                    with open(gjf_file_path, 'w') as gjf_file:
                        gjf_file.write(file_content)
                    shutil.copyfile(chk_file, new_chk_file_path)
                except FileNotFoundError:
                    print(f"Error: The checkpoint file {chk_file} was not found.")
                except IOError as e:
                    print(f"IO Error: {e}")


# Main execution
if __name__ == "__main__":
    config = read_config('config.json')

    path_read = config['path_read']
    path_write = config['path_write']
    molecule_name = config['molecule_name']
    file_ext_gjf = config['file_ext_gjf']
    file_ext_log = config['file_ext_log']
    file_ext_chk = config['file_ext_chk']
    normal_modes = config['normal_modes']
    displacements = config['displacements']

    input_gjf_file = os.path.join(path_read, f'{molecule_name}_0{file_ext_gjf}')
    output_log_file = os.path.join(path_read, f'{molecule_name}_0{file_ext_log}')
    initial_chk_file = os.path.join(path_read, f'{molecule_name}_0{file_ext_chk}')

    normal_mode_coordinates, num_atoms = read_normal_mode_displacements(output_log_file)
    if normal_mode_coordinates is None or num_atoms is None:
        exit(1)

    atom_labels, x_coordinates, y_coordinates, z_coordinates, header_lines, tail_lines = read_ground_state_coordinates(
        input_gjf_file, num_atoms)
    if atom_labels is None or x_coordinates is None:
        exit(1)

    write_gjf_files(path_write, molecule_name, normal_modes, displacements, atom_labels, x_coordinates, y_coordinates,
                    z_coordinates, normal_mode_coordinates, header_lines, tail_lines, initial_chk_file)
