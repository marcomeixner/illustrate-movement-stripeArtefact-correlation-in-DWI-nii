import os
import subprocess
import numpy as np
import shutil
import matplotlib.pyplot as plt
import nibabel as nib
import csv
from scipy.fftpack import fft2, fftshift
import glob


# define function to extract movement paraemters from affine matrix
def extract_affine_components(output_path):
    """
    Extract translation and scaling components from affine matrices stored in a directory.

    Args:
        output_path (str): Path to the base directory containing affine matrices in "MAT_" subfolder.

    Returns:
        tuple: Six numpy arrays containing:
            - tx (array): Translation in X.
            - ty (array): Translation in Y.
            - tz (array): Translation in Z.
            - scale_x (array): Scaling in X.
            - scale_y (array): Scaling in Y.
            - scale_z (array): Scaling in Z.
    """
    
    
    if not os.path.isdir(output_path):
        raise FileNotFoundError(f"Affine matrices directory not found: {affine_matrices_dir}")


    # Path to the directory containing affine matrices
    affine_matrices_dir = os.path.join(output_path, "eddy.eddy_parameters")

    # Initialize lists to store vector components
    tx, ty, tz = [], [], []    
    rot_x, rot_y, rot_z = [], [], []

    affine_data = np.loadtxt(affine_matrices_dir, usecols=range(6))
    tx, ty, tz, rot_x, rot_y, rot_z = affine_data.T  # Transpose to unpack columns

    # Convert to numpy arrays
    tx = np.array(tx)
    ty = np.array(ty)
    tz = np.array(tz)
    
    rot_x = np.array(rot_x)
    rot_y = np.array(rot_y)
    rot_z = np.array(rot_z)
    
    return tx, ty, tz, rot_x, rot_y, rot_z
    


# define function to extract a score for how 'stripy' a sagittal projection of the nii inputj image is
def get_stripe_match_score(twoD_image, periodicity, bandwidth):

    # Perform Fourier Transform
    F = fft2(twoD_image)
    F_shifted = fftshift(F)
    magnitude_spectrum = np.abs(F_shifted)

    # Identify the target frequency index
    freq_index = int(twoD_image.shape[0] / periodicity)

    # Create a mask to isolate the target frequency band
    bandwidth = 5  # Allow a small range around the target frequency
    mask = np.zeros_like(F_shifted, dtype=bool)
    mask[freq_index - bandwidth:freq_index + bandwidth, :] = True

    # Compute the total energy and energy in the target frequency band
    total_energy = np.sum(magnitude_spectrum)
    target_energy = np.sum(magnitude_spectrum[mask])

    # Compute the match score
    match_score = target_energy / total_energy
    
    return match_score



###############################################################################################
# Define input and output paths here
###############################################################################################


# project specific:
input_path = "/mnt/c/Users/marco/sciebo/code/DWI_sagittal_artefact/01_github/niiData"
    # contains nii (or nii.gz) files and acq parameter files for fsl eddy (acqparams.txt, bval.bval, bvec.bvec, index.txt) 

output_path_proj = '/mnt/c/Users/marco/sciebo/code/DWI_sagittal_artefact/01_github/out/github'
fsl_output_path = "/mnt/c/Users/marco/sciebo/code/DWI_sagittal_artefact/01_github/fsl_out/"  # Path for output 4D file

inputFormat = '*.nii'
periodicity = 2

# make a list of all nii files found in the input dir
nii_files = glob.glob(os.path.join(input_path, inputFormat))
print('all nii files found:', len(nii_files))

if 1 == 1: # option to check which txt are already in the output and only run through those nii-files that are not
    # isolate nii filenames 
    nii_filenames = [os.path.splitext(os.path.basename(path))[0] for path in nii_files]
    nii_filenames = [s.lstrip("_") for s in nii_filenames]
    # isolate txt filenames
    txt_files = glob.glob(os.path.join(output_path_proj, '*.txt'))
    txt_filenames = [os.path.splitext(os.path.basename(path))[0] for path in txt_files]
    # get the delta of nii and txt filenames
    delta = list(set(nii_filenames) ^ set(txt_filenames))
    # filter for the delta (Keep only elements from nii_files that contain any entry of delta as a substring)
    nii_files = [item for item in nii_files if any(sub in item for sub in delta)]

    print('unprocessed nii files found:', len(nii_files))



# loop over all nii files:
###################################
cnt = 1
for nii_file_path in nii_files:


    # isolate the nii filename and remove '_':
    nii_filename = os.path.basename(nii_file_path)  # Extract only the filename
    nii_filename = nii_filename.lstrip('_')
    nii_filenameNoExt, extension = os.path.splitext(nii_filename)
    if inputFormat == '*.nii.gz':
        nii_filenameNoExt, extension = os.path.splitext(nii_filenameNoExt)


    output_path = os.path.join(output_path_proj, nii_filenameNoExt + '.txt')

    # print out of each run
    print('###########################################################################')
    print(f"Run: {cnt}")
    print(f"Periodicity: {periodicity}")    
    print(f"Processing file path: {nii_file_path}")
    print(f"Processing file name: {nii_filename}")
    print(f"out path: {output_path}")
    print(f"fsl_output_path: {fsl_output_path}")
    print('###########################################################################')

    cnt += 1

    # empty fsl_output_path
    try:
        # Loop through all items in the folder
        for item in os.listdir(fsl_output_path):
            item_path = os.path.join(fsl_output_path, item)

            # Check if the item is a directory
            if os.path.isdir(item_path):
                shutil.rmtree(item_path)  # Remove directory and its contents
                print(f"Deleted directory: {item_path}")
            else:
                os.remove(item_path)  # Remove file
                print(f"Deleted file: {item_path}")

    except Exception as e:
        print(f"Error while emptying folder: {e}")



    # use fsl bet to make a mask, which is required for fsl eddy:
    ################################################################
    mask_path = os.path.join(fsl_output_path, "mask.nii")
    bet_command = ["bet", nii_file_path, mask_path, "-m"]
    subprocess.run(bet_command, check=True)
    
    
    
    
    # Set the FSL eddy command (expected to be a list of strings; the first item is the command, the following items are its arguments)
    #############################
  
    index_path = os.path.join(input_path, "acq_par", "index.txt") # enter col of ones according to num vols
    acqp_path = os.path.join(input_path, "acq_par", "acqparams.txt")   # enter total read out time
    bvecs_path = os.path.join(input_path, "acq_par", "bvec.bvec") # rename to bvec
    bvals_path = os.path.join(input_path, "acq_par", "bval.bval")   # rename to bval
    eddy_path = os.path.join(fsl_output_path, "eddy")

    eddy_command = [
        "eddy",
        f"--imain={nii_file_path}",
        f"--mask={mask_path}",
        f"--index={index_path}",
        f"--acqp={acqp_path}",
        f"--bvecs={bvecs_path}",
        f"--bvals={bvals_path}",
        f"--out={eddy_path}"
    ]
    
    subprocess.run(eddy_command, check=True)
    print(f"Motion correction completed. Output saved to: {fsl_output_path}")

    # extract affine components
    tra_x, tra_y, tra_z, rot_x, rot_y, rot_z = extract_affine_components(fsl_output_path)

    # calculate the change in movement
    diff_tra_z = [0] 
    for i in range(1, len(tra_z)):
        #print(i)    
        #print(tra_z[i])
        if tra_z[i - 1] < tra_z[i]:
            diff_tra_z.append(abs(tra_z[i - 1] - tra_z[i]))
        else:
            diff_tra_z.append(-abs(tra_z[i - 1] - tra_z[i]))
            

    # get match score for nifti 
    nifti_image = nib.load(nii_file_path)   # returns a NIfTI image object that contains both the image data and associated metadata
    nifti_data = nifti_image.get_fdata() # retrieves the image data stored in the NIfTI file and converts it into a NumPy array
    nifti_data_sagProj = np.sum(nifti_data, axis=0)

    nii_proj_match_score = []

    bandwidth = 5

    for volIdx in range(nifti_data.shape[3]):
        nii_proj_match_score.append(get_stripe_match_score(np.rot90(nifti_data_sagProj[:, :, volIdx]), periodicity, bandwidth))



    # save result to outpath:
    #################################
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')  # Use tab as delimiter for a neat table
        #writer.writerow(['Column 1', 'Column 2'])  # Writing header (optional)
        for v1, v2 in zip(diff_tra_z, nii_proj_match_score):
            writer.writerow([v1, v2])  # Writing data rows


    