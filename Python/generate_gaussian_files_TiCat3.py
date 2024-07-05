import re

# Define paths and filenames (modify as needed)
source_path = "C:/Users/Thomas/Desktop/tic"
destination_path = "C:/Users/Thomas/Desktop/tic"
source_filename = "TiCat3_1"
destination_filename_prefix = "TiCat3_"
file_extension = ".gjf"

# Read source file contents
with open(f"{source_path}/{source_filename}{file_extension}", "r") as source_file:
    file_contents = source_file.read()

# Loop through states (2 to 20, inclusive)
for state in range(2, 21):
    # Modify contents with state number
    updated_contents = re.sub(r"%chk=TiCat3_1.chk", f"%chk=TiCat3_{state}.chk", file_contents)
    updated_contents = re.sub(r"nstates=20,root=1\)", f"nstates=20,root={state})", updated_contents)

    # Construct destination filename
    destination_filename = f"{destination_path}/{destination_filename_prefix}{state}{file_extension}"

    # Write modified contents to destination file
    with open(destination_filename, "w") as destination_file:
        destination_file.write(updated_contents)

