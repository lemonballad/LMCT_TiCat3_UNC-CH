import re

# Define paths and filenames (modify as needed)
source_path = "C:/Users/Thomas/Desktop"  # Example path
old_path = source_path  # Same as source directory
new_path = source_path  # Example path

source_filename = "TiCat3"
old_filename = f"{source_filename}_old"
new_filename = f"{source_filename}_new"
file_extension = ".log"  # File extension


# Function to construct full file paths
def get_full_path(path, filename):
    if path:
        return f"{path}/{filename}{file_extension}"
    else:
        return f"{filename}{file_extension}"


# Get full file paths
source_file = get_full_path(source_path, source_filename)  # Example using function
old_file = get_full_path(old_path, old_filename)
new_file = get_full_path(new_path, new_filename)

# Check if filenames are empty
skip_correction = any(not name for name in [source_filename, file_extension])

# Copy and correct contents
if not skip_correction:
    with open(source_file, "r") as source_file:
        file_contents_old = source_file.read()

    file_contents_new = re.sub(r"Atom AN", "Atom AN", file_contents_old)

    with open(old_file, "w") as old_file:
        old_file.write(file_contents_old)

    with open(new_file, "w") as new_file:
        new_file.write(file_contents_new)

