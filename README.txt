
This software package titled: DTHB3D_Reg, performs 3D image registration for a pair of images, either synthetic or medical, which are of the same size. The user should put the images in the MedicalImages3D folder while running the code. 


The code is generated on a Windows 64-bit machine with MATLAB 2016b version. There are some functions that use parallel computation to accelerate computation. This is done by using ‘parfor’ loops. Also, we have use mex functions which also accelerate the code. 

General instructions to run the code successfully:
1. I1 = Moving image and I2 = Fixed image
2. The default number of maximum refinement levels are set as 3. Code will work for refinement level from 1 to 4. 
3. To run the files, go to the run folder and run the appropriate mainfunction_*.m files for each of the example listed below.
4. Before running the mainfunction_*.m files, first run the compile_mex.m file in the run/ folder. This function compiles all the MEX functions to run in parallel in OpenMP framework. 



For different examples, we have used different main files for running with the respective settings:

1. Sphere to torus image example (Fig. 2 in article)
run file: mainfunction_synthetic_img.m function is used to run this example. In this function the .mat files should be correctly loaded in lines 35 onwards, i.e. sphere_200_255.mat and torus_200_255.mat files. The image files are in synthetic_images folder. User should make sure I1 = sphere_200_255 and I2 = torus_200_255. The output file containing the resulting output of the accuracy, number of control points, etc is written in output_file_sphere_torus.txt file. Modify filename of the output file as per the example.

The error tolerance value is set as tol = 10e-04 and the maximum iterations are set as 50. Users can set it accordingly on lines 30 onwards. The rest of the parameters are set in the file setparameters_synthetic.m file in the setparameter_files folders. The function call is on line 20 in the main function file.

To set the parameters, modify the setparameters_synthetic.m function. 

2. Sphere to Sun example (Fig. 3 in article):
run file: mainfunction_synthetic_img.m function is used to run this example. The image files are in synthetic_images folder. In this function the .mat files should be correctly loaded in lines 35 onwards, i.e. Sphere.mat and sun_like.mat files. User should make sure I1 = Sphere and I2 = sun_like. The output file containing the resulting output of the accuracy, number of control points, etc is written in output_file_sphere_sun.txt file. Modify filename as per the example. To set the parameters, modify the setparameters_synthetic.m function.

3. Brain MRI taken from Brainweb website (Fig. 5 in article)
run file: mainfunction_medical_img.m function is used to run this example. Download the files from the website (http://brainweb.bic.mni.mcgill.ca/brainweb/anatomic_normal_20.html). Save the files in medical_images folder. The output file containing the resulting output of the accuracy, number of control points, etc is written in output_file_medical_img.txt file. Modify filename as per the example. To set the parameters, modify the setparameters.m function.

4. Registration of pre-operative and post-operative MRI (Fig. 9)
run file: mainfunction_medical_img_prepost.m function is used to run this example.  Download the dataset of the resection surgery from the website: (http://www.spl.harvard.edu/publications/item/view/1915). The zip file ResectionDataset2010.zip is to be extracted. We have the pre-operative image volume (PreOp.img and PreOp.hdr) and post-operative (PreOp.img and PreOp.hdr). Store all the four files ie. PreOp.hdr, PreOp.img, PostOp.hdr and PostOp.img files in medical_images folder. There will a guy window opened and under the drop down menu of the ReadData3D window select under ‘Image file format’ ad ‘HDR/IMG Analyze (*hdr)’ and click Load. Similarly do that for PostOp.img file. The read3DImage_prepost.m file stores the image marines and makes the images of equal size. The images are then used as the input for the image registration algorithm. The output file containing the resulting output of the accuracy, number of control points, etc is written in output_file_prepost.txt file on line 26 of the main function file. Modify filename as per the example. To set the parameters, modify the setparameters_prepost.m function on line 17.


5. Registration of liver CT images (Fig. 8)
run file: mainfunction_medical_img_liverCT.m function is used to run this example.  The image files are to be downloaded from (http://www.insight-journal.org/midas/collection/view/38.) Click on Original Datasets and download Patient01.mha and Patient04.mha files and save them in the medical_images folder. The read3Dimage_liverCT.m file reads the .mha files and resizes the images to make the sizes equal. Rigid alignment of the images is carried out using imregister() function  in MATLAB with ‘rigid’ option. The output image matrices are then registered in the non-rigid framework in our code. The output file containing the resulting output of the accuracy, number of control points, etc is written in output_file_liver.txt file. Modify filename as per the example. To set the parameters, modify the setparameters_liver.m function given on line 17.

6. Registration of brain MR T1 set of images (Fig. 7)
run file: mainfunction_medical_img_brainT1.m function is used to run this example. The images of two subjects 105 and 109 were taken from (http://www.insight-journal.org/rire/download_data.php). The mr_T1 files were downloaded and converted to readable format as the input to the code. The files to be used are patient_105_mr_T1.mha and patient_109_mr_T1.mha. The image files are already present in the converted .mha format for easy use. There is no need to download the files from the website. The output file containing the resulting output of the accuracy, number of control points, etc is written in output_file_brainMR_T1.txt file. Modify filename as per the example. To set the parameters, modify the setparameters_mr_T1.m function given on line 17.