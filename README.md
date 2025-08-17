# Kikhia et al 2025
This repository hosts the code used in the manuscript: "Multicolor fate mapping of microglia reveals polyclonal expansion, heterogeneity, and cell-cell interactions after ischemic stroke"

# Content
* Each subdirectory contains the code files necessary to run the analysis indicated in the subdirectory's name.

# Data
* To run this code, please download the associated data available here: https://osf.io/h9fjy/, which contains the following subdirectories:
	* "Data" contains the following: 
		* Cell location for one animal per group for DBSCAN and Monte Carlo Simulation demos.
		* 4 raw z-stacks with their U-Net segmentation results for the morphology analysis. 
		* The spreadsheet for the Ki-67 analysis.
	* "Microglia_Morphology" contains the pretrained U-Net model.
	* "Monte_Carlo" contains the intermediary output of the analysis to save time.  
* To ensure that the code finds the data files, the data files need to be stored in the right paths. "FILE_TREE.txt" shows the correct locations of the files.  

# System requirements
* For cluster analysis and Monte Carlo simulation:
	* MATLAB R2022a with the following toolboxes:
		* Statistics and Machine Learning Toolbox
		* Antenna Toolbox™
		* Parallel Computing Toolbox
* For morphological analysis: 
	* Python 3.9 with the following libraries:
		* Jupyter Notebook 7.0.6   
		* TensorFlow-GPU 2.11.0
		* Keras 2.11.0
		* NumPy 1.24.2
		* patchify 0.2.3 
		* glob2 0.7
		* czifile 2019.7.2 
		* tifffile 2023.2.28
		* Matplotlib 3.7.0
		* Pandas 2.1.4
		* seaborn 0.13.2
		* pingouin 0.5.4
		* scikit_posthocs 0.8.1 
		* ipywidgets 8.1.1
		* scikit-learn 1.6.1
		* umap-learn 0.5.7   
	* Fiji (ImageJ 1.54f) with the following plugins:
		* MorphoLibJ. Website: https://imagej.net/plugins/morpholibj. Plugin name in the update list: IJPB-plugins. 
		* 3D ImageJ Suite. Website: https://imagej.net/plugins/3d-imagej-suite/.

# Hardware
* The code files were run on a standard laptop with the following configurations:
	* Processor: AMD Ryzen 7 5825U with Radeon Graphics   2.00 GHz
	* Installed RAM: 16.0 GB (15.4 GB usable)
* Training and application of the U-Net model required GPU acceleration. They were performed on an HPC cluster.

# Installation guide
* Download the source code with the demo data. 
* Install MATLAB with the dependencies mentioned above.
* Create a virtual environment with Python and the libraries mentioned above
* Install Fiji and include the dependencies in the update list.

# Installation time: 
* Downloading the code with the demo data requires several minutes, depending on the internet speed.
* Installing all necessary software and packages might require a couple of hours. However, the code can be run immediately if all required software is available.

# Instructions for demo analysis
After downloading the required dependencies, the demo data allows running the following analysis: 
* DBSCAN (runtime =~ 20 seconds):
	* Open the code file and run the code.
	* Output: a table containing the numerical values underlying the cluster analysis presented in Fig.2 e, f. The demo data contains one animal per group.
* Cluster analysis with Ki-67 (runtime =~ 5 seconds):
	* Open the code file "Ki67_analysis" and run the code.
* Monte Carlo (runtime =~ 2 hours)
	* Open the code file Monte_Carlo_a, and run the code.
	* Then run the codes Monte_Carlo_b, and Monte_Carlo_c sequentially. 
	* Output: a one-animal-per-group version of Fig. 2c. Intermediary outputs from a previous run are saved in the directory "Monte_Carlo/Output/" in the data repo. Using the intermediary outputs, the code can be run in a few minutes. See instructions in the code files.  
* Morphological analysis (runtime =~ 20 minutes):
	* To segment RAW czi files using the pretrained U-Net model:
		* Open the Juypter Notebook "unet_model_application_smoothed_edges.ipynb". Make sure the input directory to the raw data is correct "../Data/Morphology_analysis/RAW/". Run the notebook. 
		* The code will create tif files with binary images as a result of the segmentation. The same results are already saved from a previous run in the data repo at "Data/Morphology_analysis/Unet_output/".
	* To segment single cells and extract their morphological features similar to Fig. 5a:
		* Open Fiji and open the ImageJ macro code "morphology_analysis_after_unet.ijm"
		* Make sure to change the input directories for both RAW images and U-Net output images. 
		* Run the code.
		* To produce the exact same results in Fig. 5a, provide the code with the IDs of excluded nuclei when prompted. The information is available in the file "excluded_nuclei_demo_data.txt" 
		* Output: the images shown in Fig. 5a, and CSV files for further statistical analysis.                        

# Licenses
* The codes in this repository are licensed under the BSD-style license found in the LICENSE file in the root directory of this source tree.
* The code smooth_tiles_predictions.py is licensed under an MIT license. See more information in the code.
* The codes for training the U-Net model and applying it contain elements adapted from other authors and published under MIT license. See more information in the code files. 

# Citation
This repository accompanies a paper currently in preparation. See the repository’s latest version for updated references.

# Author
Majed Kikhia, Charité – Universitätsmedizin Berlin, Department of Experimental Neurology  
Codes adapted from other authors are indicated within the code files.
