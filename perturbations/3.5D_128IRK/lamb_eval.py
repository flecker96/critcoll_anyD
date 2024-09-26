import os
import shutil
import subprocess

# Base value
lambguess = 2.59

for lamb in (lambguess,0.0):

    #Output files of shoot_main.exe to be copied
    output_files=(f"lamb.dat",f"doutdin.junk",f"status.junk",f"detoflamb.junk")

    # Name of the .dat file to be created
    lamb_file = f"lamb.dat"

    # Name of the new directory to be created
    new_directory = f"lambval_{lamb}"

    # Write the content to the .dat file
    with open(lamb_file, 'w') as f:
        f.write(str(lamb))

    print(f"File '{lamb_file}' has been created with: lamb = {lamb}\t")

    # Run the shoot_inner.exe file
    evolve = subprocess.Popen(["./shoot_inner.exe"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Iterate over the output lines in real time
    for line in iter(evolve.stdout.readline, ''):
        print(line, end='')

    # Wait for the shooting to complete
    evolve.stdout.close()
    evolve.wait()

    # Create the new directory
    os.makedirs(new_directory, exist_ok=True)

    # Check if the directory was created successfully
    if os.path.exists(new_directory) == False:
        print(f"Failed to create directory '{new_directory}'.")
        exit(1)

    # Copy the output files to the new directory
    for f in output_files:
        shutil.copy(f, new_directory)
        # Check if the copy operation was successful
        if os.path.exists(os.path.join(new_directory, f)) == False:
            print(f"Failed to copy file '{f}' to '{new_directory}'.")
            exit(1)

    for f in output_files:
        os.remove(f)

