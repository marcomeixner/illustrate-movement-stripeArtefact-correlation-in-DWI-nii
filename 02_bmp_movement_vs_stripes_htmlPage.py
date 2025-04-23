
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import nibabel as nib
# import cv2
import shutil
from pathlib import Path

import imageio
from PIL import Image, ImageDraw, ImageFont

# define a list of thresholds, which defines which vols are recognized as stripy (for most outer loop)
std_factor_list = [1.75] # 1.5, 1.75, 2
std_factor_list = [float(x) for x in std_factor_list]


# list of all projects (for next loop)
projects = ["github_proj"]   # comma seperated list of projects
niiFormatList = ['.nii']   # comma seperated list:  '.nii','.nii.gz'
nii_path_list = ['/mnt/c/Users/marco/sciebo/code/DWI_sagittal_artefact/01_github/niiData']  # comma seperated list


# example of three projects
'''
# specific projects
projects = ["proj_1", "proj_2", "proj_3"]
niiFormatList = ['.nii', '.nii', '.nii.gz']
nii_path_list = ['/path/to/nii/files/of/proj/1',
                 '/path/to/nii/files/of/proj/2',
                 '/path/to/nii/files/of/proj/3']
'''

# decide if you want to delete all existing html-output (1 is yes)
del_HTML_out = 0

# where are the z-movement vs stripe txts:
txt_parent_path = '/mnt/c/Users/marco/sciebo/code/DWI_sagittal_artefact/01_github/out'

# where to put the html output:
html_output_parent = "/mnt/c/Users/marco/sciebo/code/DWI_sagittal_artefact/01_github/"
# html_output_dir = "/mnt/c/Users/marco/sciebo/code/DWI_sagittal_artefact/html_out"


# decide if you want to generate a plot for the website and stuff:
enablePlot = 1


def generate_html(directory, projTitle, std_factor, proj_avg_corr, proj_outlier_Pcnt, proj_num_outlier, num_datasets, proj_total_vols_num, unq_bvals):
    folder = Path(directory)
    images = sorted(folder.glob("*.png"))  # Get all PNG files
    html_filename = projTitle + ".html"
    output_path = os.path.join(folder, html_filename)
    
    unq_bvals_str = ', '.join(map(str, unq_bvals))    

    # make list of outlier datasets to add to html
    outlier_data = [str(img) for img in images]
    outlier_data = outlier_data[0::2]
    outlier_data = [os.path.splitext(os.path.basename(path))[0] for path in outlier_data]
    outlier_data = [item.rstrip('_imshow') for item in outlier_data]
    num_outlier_ds = len(outlier_data)
    #print("outlier_data: ", outlier_data)

    
    # Default text positions
    text_positions = {}
        
    html_content = f"""<!DOCTYPE html>
    <html lang='en'>
    <head>
        <meta charset='UTF-8'>
        <meta name='viewport' content='width=device-width, initial-scale=1.0'>
        <title>{projTitle}</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                text-align: center;
                background-color: black;  /* Set the background to black */
                color: white;             /* Set the text color to white */
            }}
            .gallery {{
                display: flex;
                flex-direction: column;  /* Ensure content stays vertical */
                align-items: center;
                gap: 20px;
            }}
            .gif-gallery {{
                display: flex;
                justify-content: center; /* Center the GIFs horizontally */
                gap: 20px;               /* Space between GIFs */
            }}
            .image-container {{
                position: relative;
                display: block;
                margin-bottom: 20px;
            }}
            img {{
                max-width: 250px;    /* Adjust the size of GIFs */
                height: auto;        /* Maintain aspect ratio */
                display: block;
            }}
            .text-overlay {{
                position: absolute;
                background: rgba(0, 0, 0, 0.6);
                color: white;
                padding: 5px;
                border-radius: 5px;
            }}
            .image-label {{
                font-size: 18px;
                margin-bottom: 10px;
            }}
            .hidden {{
                display: none;
            }}
        </style>
        <script>
            function checkPassword() {{
                var password = prompt('Please enter the password:');
                if (password === 'password') {{
                    document.getElementById('protected-content').classList.remove('hidden');
                }} else {{
                    alert('Incorrect password.');
                    document.body.innerHTML = 'Access denied.';
                }}
            }}
            window.onload = checkPassword;
        </script>
    </head>
    <body>
        <h1>{projTitle}</h1>
        <div id="protected-content" class="hidden">
            <div class='gallery'>
                <p style='font-size: 20px;'><strong>In total {num_datasets} datasets were checked for a stripe pattern. For a threshold of {std_factor} (std_factor) stripes were found in {num_outlier_ds} datasets.</strong></p>            
                <p style='font-size: 20px;'>std_factor: {std_factor}</p>
                <p style='font-size: 20px;'>unique b-values: {unq_bvals_str}</p>
                <p style='font-size: 20px;'>total_vols_num: {proj_total_vols_num}</p>
                <p style='font-size: 20px;'>num_outlier: {proj_num_outlier}</p>
                <p style='font-size: 20px;'>outlier_Pcnt: {proj_outlier_Pcnt}%</p>
                <p style='font-size: 20px;'>avg correlation: {proj_avg_corr}</p>
                <p style='font-size: 20px;'><strong>datasets with outliers detected:</strong><br><br>{'<br>'.join(outlier_data)}</p>
                """

    # loop over all images found in the folder:
    for index, img in enumerate(images):
        img_name = img.name
        text = text_positions.get(img_name, "")

        img_name_noExt = Path(img_name).stem
        img_name_nii = Path(img_name_noExt).with_suffix(".nii")

        if index % 2 == 0:  # Add text above every second image (starting with the first)
            html_content += f"<p style='font-size: 35px;'><strong>File: {img_name_nii}</strong></p>"

            # ADD GIFS
            ##############

            # Wrap the GIFs in a container that will arrange them horizontally
            html_content += "<div class='gif-gallery'>"

            # add all vol gif:
            gif_name_allVol = img_name_noExt + ".gif"
            gif_path_allVol = "bval_gifs/" + gif_name_allVol
            gif_path_allVol = gif_path_allVol.replace("imshow", "allVols")
            label_text_allVol = "all Volumes"  # Your label as a variable

            html_content += f"""
            <div class='image-container'>
                <span class="image-label"><strong>{label_text_allVol}</strong></span>
                <img src='{gif_path_allVol}' alt='GIF Image' width="250">
                {f"<div class='text-overlay'>{text}</div>" if text else ""}
            </div>
            """
        
            # add bval gifs:
            #################
            
            # max show 4 bvals to the page does not get too wide
            if len(unq_bvals) > 4:
                unq_bvals = unq_bvals[:4]  # Keep only the first 4 elements
            
            for bvalIdx in range(len(unq_bvals)):
                gif_name = img_name_noExt + str(unq_bvals[bvalIdx]) + ".gif"
                gif_path = "bval_gifs/" + gif_name
                gif_path = gif_path.replace("imshow", "bval")
                label_text = "bval: " + str(unq_bvals[bvalIdx])  # Your label as a variable

                html_content += f"""
                <div class='image-container'>
                    <span class="image-label"><strong>{label_text}</strong></span>
                    <img src='{gif_path}' alt='GIF Image' width="250">
                    {f"<div class='text-overlay'>{text}</div>" if text else ""}
                </div>
                """

            html_content += "</div>"  # Close the gif-gallery div


        # if index % 2 != 0:   # option to only plot the graphs, not the nii slices

        html_content += f"""
        <div class='image-container'>
            <img src='{img_name}' alt='{img_name}' style="max-width: 1200px; height: auto;">
            {f"<div class='text-overlay'>{text}</div>" if text else ""}
        </div>
        """




    html_content += """
            </div>
        </div>
    </body>
    </html>
    """





    
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    
    print(f"HTML gallery saved as {output_path}")


for std_factor in std_factor_list:   # loop over std values

    html_fol_name = f"html_out_{std_factor}"
    html_fol_path = os.path.join(html_output_parent, html_fol_name)


    print("std_factor_list: ", std_factor_list)
    print("html_fol_path: ", html_fol_path)


    # Check if the folder exists
    if os.path.exists(html_fol_path):
        # If the folder is not empty, use shutil.rmtree() to delete the folder and all its contents
        shutil.rmtree(html_fol_path)
        print(f"The folder '{html_fol_path}' has been deleted.")
    else:
        pass


    # make parent html out folder
    if del_HTML_out == 1: # delete existing even if there is one
        os.mkdir(html_fol_path)     # always makes a new folder
    else: # only makes a new folder if there is none
        os.makedirs(html_fol_path, exist_ok=True) 


    # copy index.html
    #index_path = os.path.join(html_output_parent, "index.html")
    #shutil.copy(index_path, html_fol_path)
    # copy compareDWIseq.xlsx
    #compareDWIseq_path = os.path.join(html_output_parent, "compareDWIseq.xlsx")
    #shutil.copy(compareDWIseq_path, html_fol_path)       


# html_output_parent = "/mnt/c/Users/marco/sciebo/code/DWI_sagittal_artefact"
    for projName, nii_path, niiFormat in zip(projects, nii_path_list, niiFormatList):   # loop over projects

        # read in bvals and order:
        bval_path = os.path.join(html_output_parent, "niiData", "acq_par", "bval.bval")
        with open(bval_path, 'r') as file:
            for line in file: # Loop through each line in the file
                bvals = line.strip().split() # Strip any leading/trailing whitespace and split the line by spaces
                bvals = [int(value) for value in bvals]  # Convert each value from string to integer

        unq_bvals = list(set(bvals))
        unq_bvals.sort()
        
        bval_dict = {} # Initialize an empty dictionary to store indices for each value
        for index, value in enumerate(bvals): # Loop through the vector and store indices for each unique value
            if value not in bval_dict:
                bval_dict[value] = []  # Initialize an empty list for the value
            bval_dict[value].append(index)
        for value in bval_dict: # Iterate through each unique number
            bval_for_value = bval_dict.get(value, []) 

        print('##################################')
        print(projName, ": ", "std_factor ", std_factor)
        print("unq_bvals: ", unq_bvals)
        print('##################################')



        html_path = os.path.join(html_fol_path, projName)

        # make project html out folder
        os.mkdir(html_path)
     

        
        # delete existing folder 
        ##############################
        if enablePlot:
            try:
                if os.path.exists(html_path) and os.path.isdir(html_path):
                    shutil.rmtree(html_path)
            except:
                pass  # Silently handle any errors without printing
            
            # make new html output dir
            os.makedirs(html_path, exist_ok=True)




        inputDir = os.path.join(txt_parent_path, projName)


        txt_files = glob.glob(os.path.join(inputDir, '*.txt'))
        # print("txt_files: ", txt_files)
        
        num_datasets = len(txt_files)  # Count the number of txt files


        corr_vec = []
        corr_vec_abs = []
        
        fullOutlierList = []
        numberOfVols = 0

        cnt = 0

        for txt_path in txt_files:  # loop txt files (eg datasets)

            data = np.loadtxt(txt_path, delimiter="\t")


            z_movement = data[:,0]
            match_score = data[:,1]

            numberOfVols = numberOfVols + len(z_movement)


            
    
            cnt += 1

    
            # get deviation vector:
            scores = [float(score) for score in match_score] # Convert to floats
            std_factor = float(std_factor)
            mean_value = np.mean(scores)
            std_value = np.std(scores)
            threshold = std_factor * std_value
            deviation_vector = [score - mean_value for score in scores]
            deviation_vector_abs = np.abs(deviation_vector)    
            
            
            corr_vec.append(np.corrcoef(z_movement, deviation_vector)[0, 1])
            
            curr_corr_abs = np.corrcoef(abs(z_movement), deviation_vector_abs)[0, 1]
            corr_vec_abs.append(curr_corr_abs)        


            # get match score outliers:    
            match_score_mean = np.mean(match_score)
            match_score_std = np.std(match_score)
            # Identify outliers as those beyond the threshold of mean ± (std_factor * std_value)
            match_score_outliers = [
                i for i, score in enumerate(match_score) if abs(score - match_score_mean) > std_factor * match_score_std
                ]


            # remove first vol from outliers if present
            if 0 in match_score_outliers:
                match_score_outliers = [x for x in match_score_outliers if x != 0]  # Remove 0 only if present
                
        
            # check if there are any match_score_outliers:
            if match_score_outliers:


                # get z-movement outliers:
                z_movement_mean = np.mean(z_movement)
                z_movement_std = np.std(z_movement)
                # Identify outliers as those beyond the threshold of mean ± (std_factor * std_value)
                z_movement_outliers = [
                    i for i, score in enumerate(z_movement) if abs(score - z_movement_mean) > std_factor * z_movement_std
                    ]



                # full list of outliers:
                if enablePlot:
                    print(os.path.basename(txt_path), ": ", match_score_outliers)  


                # create the nii path:
                txt_name = os.path.basename(txt_path)      
                txt_name_only = os.path.splitext(txt_name)[0]
                
                nii_path_file = os.path.join(nii_path, txt_name_only) + niiFormat  

                if enablePlot:         
                    print("nii_path_file: ", nii_path_file)

                num_match_outlier = len(match_score_outliers)

                numSagImg = 5
                    
                
                vert_merged_image = None


                # create a gif files according to bvals: 
                ########################################################################################
                nii_data = nib.load(nii_path_file)
                nii_4d = nii_data.get_fdata()
                niiShape = nii_4d.shape

                # 1st: create for the whole volume:
                ###################################
                images_all = []
                for volum in range(niiShape[3]): # loop over volumes of current bval to make the gif
                    img = nii_4d[int(np.round(niiShape[0]/2)),:,:,volum]
                    img = np.rot90(img) # # rotate image
                    img_normalized = (img - np.min(img)) / (np.max(img) - np.min(img)) * 255 # Normalize the image to 0-255
                    img_normalized = img_normalized.astype(np.uint8) # Convert to uint8
                    image_pil = Image.fromarray(img_normalized)   # Convert the np array
                    draw = ImageDraw.Draw(image_pil)  # add label to upper right corner
                    label = f"{volum + 1}, {bvals[volum-1]}"  # DEFFINE LABEL HERE
                    font = ImageFont.load_default() # Load the default font
                    text_bbox = draw.textbbox((0, 0), label, font=font)  # Get bounding box coordinates
                    text_width = text_bbox[2] - text_bbox[0]  # Width of the text
                    text_height = text_bbox[3] - text_bbox[1]  # Height of the text
                    text_position = (3, 3)  # Position at bottom-right
                    draw.text(text_position, label, font=font, fill="white")
                    images_all.append(image_pil) # appen to images_all

                # Save images_all as a GIF                
                plot_gif_filename= txt_name_only + "_allVols" + ".gif"
                plot_gif_path = os.path.join(html_path, "bval_gifs")
                os.makedirs(plot_gif_path, exist_ok=True)
                plot_gif_fullpath = os.path.join(plot_gif_path, plot_gif_filename)
                imageio.mimsave(plot_gif_fullpath, images_all, duration=0.28, loop=0)  # duration is the time per frame in seconds

                # 2nd: create a gif for each bval:
                ##################################
                for bvalIdx in range(len(unq_bvals)): # loop over bvals
                    indices_for_value = bval_dict.get(unq_bvals[bvalIdx], [])
                    images = []
                    for volum in indices_for_value: # loop over volumes of current bval to make the gif
                        img = nii_4d[int(np.round(niiShape[0]/2)),:,:,volum]
                        img = np.rot90(img) # # rotate image
                        img_normalized = (img - np.min(img)) / (np.max(img) - np.min(img)) * 255 # Normalize the image to 0-255
                        img_normalized = img_normalized.astype(np.uint8) # Convert to uint8
                        image_pil = Image.fromarray(img_normalized)   # Convert the np array
                        draw = ImageDraw.Draw(image_pil)  # add label to upper right corner
                        label = f"{volum + 1}"  # DEFFINE LABEL HERE
                        font = ImageFont.load_default() # Load the default font
                        text_bbox = draw.textbbox((0, 0), label, font=font)  # Get bounding box coordinates
                        text_width = text_bbox[2] - text_bbox[0]  # Width of the text
                        text_height = text_bbox[3] - text_bbox[1]  # Height of the text
                        text_position = (3, 3)  # Position at bottom-right
                        draw.text(text_position, label, font=font, fill="white")
                        images.append(image_pil) # appen to images
                    # Save images as a GIF                
                    plot_gif_filename= txt_name_only + "_bval" + str(unq_bvals[bvalIdx]) + ".gif"
                    plot_gif_path = os.path.join(html_path, "bval_gifs")
                    os.makedirs(plot_gif_path, exist_ok=True)
                    plot_gif_fullpath = os.path.join(plot_gif_path, plot_gif_filename)
                    imageio.mimsave(plot_gif_fullpath, images, duration=0.36, loop=0)  # duration is the time per frame in seconds


                
                if enablePlot:   # plotting on / off

                
                    for volIdx in match_score_outliers: # loop over outlier volumes
                        
                        if enablePlot:
                            print('curr vol idx: ', volIdx) 
                        
                        
                        currVol = nii_data.slicer[:, :, :, volIdx].get_fdata()
                        # (sag, cor, axial)
                        #niiShape = currVol.shape          
                        #print("niiShape[0]: ", niiShape[0])
                        #print("niiShape[1]: ", niiShape[1])
                        #print("niiShape[2]: ", niiShape[2])
                        
                        
                        hor_merged_image = np.hstack([np.rot90(currVol[int(np.round(niiShape[0]/2)-np.round(numSagImg/2)+i),:,:], k=1) for i in range(numSagImg)])
                        
                        

                        # vertically merge images
                        if vert_merged_image is None:
                            vert_merged_image = hor_merged_image  # First image initializes the stack
                        else:
                            vert_merged_image = np.vstack((vert_merged_image, hor_merged_image))  # Stack new image below


                    # plot the stripe volumes:
                    plt.figure(1) 
                    plt.imshow(vert_merged_image, cmap='gray')
                    plt.axis('off')  # Hide axes
                    
                    # add labels to each vol row
                    x_pos = 5
                    y_start = niiShape[2]/2
                    step_width = niiShape[2]
                    y_positions = [y_start + i * step_width for i in range(num_match_outlier)]
                    
                    for i, y_pos in enumerate(y_positions):
                        plt.text(x_pos, y_pos, f"Vol {match_score_outliers[i]}", fontsize=(15/(num_match_outlier ** (1/3))), color="red")
                                #bbox=dict(facecol

                    # save plot one to html folder 
                    plot_one_filename= txt_name_only + "_imshow" + ".png"
                    plot_one_path = os.path.join(html_path, plot_one_filename)
                    plt.savefig(plot_one_path, bbox_inches='tight', pad_inches=0, dpi=600)
                    plt.close()
                    
                    
                    
                    
                    # subplot: movement in z-direction and matchscore outliers
                    #################################################################
                    plt.figure(2)
                    plt.figure(figsize=(20, 10))

                    # Set figure background color to black
                    plt.gcf().set_facecolor('black')
                    x = list(range(0, len(deviation_vector_abs)))

                    # Create the first subplot
                    #############################
                    plt.subplot(2, 1, 1)  # match_score subplot
                    plt.title(f"correlation of z_movement and stripes: {np.round(curr_corr_abs, 3)}", fontsize=20, color='white')  # Title in white
                    plt.plot(deviation_vector_abs, 'c-', label="match_score", lw=2)  # Plot the match score in blue
                    plt.scatter(
                        [x[i] for i in match_score_outliers],
                        [deviation_vector_abs[i] for i in match_score_outliers],
                        color='red',  # Red color for outliers
                        label=f'Outliers (> {std_factor} stds)',
                        zorder=5
                    )
                    plt.xlabel("vol idx", color='white', fontsize=20)  # X-axis label in white
                    plt.ylabel("match_score", color='white', fontsize=20)  
                    plt.axhline(y=threshold, color='r', linestyle='--', label="threshold")  # Red threshold line
                    plt.grid(True, color='white')  # White grid lines
                    
                    # Set white tick labels for both axes
                    plt.xticks(color='white', fontsize=15)  # X-axis tick labels in white
                    plt.yticks(color='white', fontsize=15)  # Y-axis tick labels in white                    

                    # Set the background of the subplot to black
                    plt.gca().set_facecolor('black')

                    # Create the second subplot
                    #############################
                    plt.subplot(2, 1, 2)  # z_movement subplot
                    plt.plot(z_movement, 'c-', label="z_movement [mm]", lw=2)  # Plot z_movement in blue
                    plt.xlabel("vol idx", color='white', fontsize=20)  # X-axis label in white
                    plt.ylabel("z_movement [mm]", color='white', fontsize=20)  
                    plt.grid(True, color='white')  # White grid lines
                    
                    # Set white tick labels for both axes
                    plt.xticks(color='white', fontsize=15)  # X-axis tick labels in white
                    plt.yticks(color='white', fontsize=15)  # Y-axis tick labels in white

                    # Set the background of the subplot to black
                    plt.gca().set_facecolor('black')
                

                    # save plot two to html folder 
                    plot_two_filename= txt_name_only + "_plot" + ".png"
                    plot_two_path = os.path.join(html_path, plot_two_filename)
                    plt.savefig(plot_two_path, bbox_inches='tight', pad_inches=0)
                    plt.close()        


                #plt.show()


                # make a full outlier list of the whole dataset for counting the total number of stripe vols
                fullOutlierList.extend(match_score_outliers)

                #print(fullOutlierList)


                #print("z_movement_outliers:", z_movement_outliers)
                #print("match_score_outliers:", match_score_outliers)



                        
            else:
                if enablePlot:
                    print("match_score_outliers is empty") 
                else:
                    pass

        if numberOfVols:
            
            
            proj_avg_corr = np.mean(corr_vec_abs)
            proj_outlier_Pcnt = round(100*len(fullOutlierList)/numberOfVols, 3)
            proj_num_outlier = len(fullOutlierList)
            proj_total_vols_num = numberOfVols
            
                    
            print("avg corr: ", proj_avg_corr)
            
            # print("avg corr: ", np.mean(corr_vec))
            # print("avg corr abs: ", np.mean(corr_vec_abs))        
            
            print("outlier-Pcnt: ", proj_outlier_Pcnt, "%")
            print("num outlier: ", proj_num_outlier)
            print("total number of vols: ", proj_total_vols_num)
            print('#########################################################')
            
    
            
            generate_html(html_path, projName, std_factor, proj_avg_corr, proj_outlier_Pcnt, proj_num_outlier, num_datasets, proj_total_vols_num, unq_bvals)


