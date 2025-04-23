This Linux Pipeline checks DWI data for a correlation between movement in z-direction (via fsl eddy) and stripe artefacts than can occur for different reasons, including movement.

The first step is to install fsl (see fsl_eddy_notes.txt).

Then the first script can be run: 01_fsl_eddy_z-movement_vs_stripeArtefact.py
In produces a txt file for each nii file. The txt file contains a line for each volume; each line shows (1) the detected movement in z-direction and (2) the 'stripyness' of the volume.

	Input:
		input_path: an Input Folder that contains nii or nii.gz files
		output_path_proj: there the correlation txt files are saved
		fsl_output_path: a Folder where the fsl eddy ouput can be stored temporarily
		inputFormat: '*.nii' or '*.nii.gz'
		periodicity: the stripe periodicity that is expected on the Images
	Output:
		a txt file for each nii file that contains a correlation table which is saved in output_path_proj


Then the second script can be run: 02_bmp_movement_vs_stripes_htmlPage.py
This script picks up the txt files, runs through the nii volumes again and creates a (Password protected) Website in html showing all nii with outliers in Terms of stripyness, including gifs for each bval.

	Input:
		std_factor_list: a list of stripyess thresholds (e.g. [1.5, 1.75, 2]) for the most outer loop
		projects: a list of Project names for the second loop
		niiFormatList: a list of Formats, one entry per Project
		nii_path_list: a list of nii paths, one for each Project
		txt_parent_path: parent path of the Output txts of the first script (containing one Folder per Project)
		html_output_parent: where the html Output Folders will go

	Output: a Folder that contains everything required for the Website, e.g.
		html_out_1.75
                └── github_proj
                    ├── bval_gifs
                    │   ├── dwi_117_AP_MaMe_20250113090131_32_allVols.gif
                    │   ├── dwi_117_AP_MaMe_20250113090131_32_bval0.gif
                    │   ├── dwi_117_AP_MaMe_20250113090131_32_bval1000.gif
                    │   ├── dwi_117_AP_MaMe_20250113090131_32_bval1800.gif
                    │   ├── dwi_117_AP_MaMe_20250113090131_32_bval2500.gif
                                    │   ├── dwi_117_AP_noMB_MaMe_20250113090131_33_allVols.gif
                    │   ├── dwi_117_AP_noMB_MaMe_20250113090131_33_bval0.gif
                    │   ├── dwi_117_AP_noMB_MaMe_20250113090131_33_bval1000.gif
                    │   ├── dwi_117_AP_noMB_MaMe_20250113090131_33_bval1800.gif
                    │   └── dwi_117_AP_noMB_MaMe_20250113090131_33_bval2500.gif
                    ├── dwi_117_AP_MaMe_20250113090131_32_imshow.png
                    ├── dwi_117_AP_MaMe_20250113090131_32_plot.png
                    ├── dwi_117_AP_noMB_MaMe_20250113090131_33_imshow.png
                    ├── dwi_117_AP_noMB_MaMe_20250113090131_33_plot.png
                    └── github_proj.html

		In case there are different Projects, the index.html file can be used to make a parent html page. Providers like netlify can be used to upload and share the Website