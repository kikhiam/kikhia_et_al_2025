// Copyright (c) 2025, Majed Kikhia
// All rights reserved.
//
// This source code is licensed under the BSD-style license found in the
// LICENSE file in the root directory of this source tree. 
//
// Author: Majed Kikhia
// February 2025
//
// Morphology analysis of microglia after segmentation with U-Net
// 
// This script is in ImageJ Macro language. It analyzes the morphology of microglial cells after performing segmentation with U-Net. 
// 
// Input: 
//		Raw images. czi files in this case.
//		U-Net segmentation output as TIF binary images.
// 		To reproduce the exat same results nuclei IDs for exclusion should be intered when prompted
//		A list of nuclei IDs are provided in "excluded_nuclei_demo_data.txt" 
// 
// The script analyzes the morphology of cells and separates touching cells using the nuclei as seeds for a 3D watershed.  
// 
// Output: 
// 		4 Tiff files: segmented and labeled cells in 3D, segmented and labelled cells in 2D with nuclei, segmented and labelled cells in 2D with skeletons overlaid, all cells in red in 2D.       
//		A CSV file containing the morphological parameters of each cell. 	
//		 
// Dependencies: 
//		MorphoLibJ. Website: https://imagej.net/plugins/morpholibj. Plugin name in the update list: IJPB-plugins.
//		3D ImageJ Suite. Website: https://imagej.net/plugins/3d-imagej-suite/. 
//
//
// clean up first
// close all images
close("*");
// empty the results table
run("Clear Results");
// configure that binary image are black in background, objects are white
setOption("BlackBackground", true);

// Setting input directories (change the directories here to match the downloaded file in your system)
Input_folder = "C:/...../MicrofettiStroke_2025/Data/Morphology_analysis/RAW";  // insert input directory for the raw images: Data/Morphology_analysis/RAW
Input_masks_folder = "C:/...../MicrofettiStroke_2025/Data/Morphology_analysis/Unet_output"; // insert input directory for the U-Net segmentation output: Data/Morphology_analysis/Unet_output

// Creating an output folder in the parent directory
Parent_folder = File.getParent(Input_folder);
Output_folder = Parent_folder + "/U-Net Morphology Analysis Output"

File.makeDirectory(Output_folder);

fileList = getFileList(Input_folder);

// Looping throw all the images
for (n=0; n < lengthOf(fileList); n++) {
	filename = fileList[n];
	Current_file = Input_folder + "/" + filename;
	
	// Opening the image	
	run("Bio-Formats Importer", "chck_for_upgrades open=[" + Current_file + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	
	// Rename the file to the file name without the whole path
	rename(filename);
	
	// Close other channeks if they exist
	getDimensions (Img_width, Img_height, Img_channels, Img_slices, Img_frames);
	getVoxelSize(width , height , depth , unit);
	
	if (Img_channels == 3) {
	run("Split Channels");
	selectWindow("C2-" + filename);
	close();
	selectWindow("C1-" + filename);	
	close();
	selectWindow("C3-" + filename);
	
	} else if (Img_channels == 2) {
	run("Split Channels");
	selectWindow("C1-" + filename);
	close();
	selectWindow("C2-" + filename);		
	}
	
	// Removing the path from the file name and informing the user about the analyzed image
	rename("DAPI_RAW");
	print("#########################################################################################");
	print("Currently processing the following image: " + filename);

	// Getting voxel size
	getVoxelSize(width , height , depth , unit);
	
	// Enhance the image of the DAPI channel
	run("Gaussian Blur...", "sigma=1 stack");
	run("Subtract Background...", "rolling=25 stack");
	run("Enhance Contrast...", "saturated=0.1 normalize process_all use");
	run("Gaussian Blur...", "sigma=1 stack");
	
	// Create a binary mask of the DAPI channel
	setAutoThreshold("Li dark stack");
	run("Convert to Mask", "method=Li background=Dark black create");
	run("Fill Holes", "stack");
	rename("DAPI_binary");
	
	// Opening U-Net mask of microglia and inform user to check that the two file names are matching
	open(Input_masks_folder + "/Segmented_" + filename + ".tif");
	
	filename_check = getTitle();
	
	print("This is the U-Net mask: " + filename_check);
	
	rename("microglia_binary");
	
	// Isolating nuclei of microglia by image calculator (AND) 
	imageCalculator("AND create stack", "microglia_binary", "DAPI_binary");
	run("Morphological Filters (3D)", "operation=Opening element=Ball x-radius=3 y-radius=3 z-radius=3");
	run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
	run("Label Size Filtering", "operation=Greater_Than size=1500");
	
	// Removing any mislabeled nuclei, cells on the border, or cells with unsatisfactory segmentation
	// To remove nothing just click OK 
	Dialog.create("Nuclei IDs to exclude");
	Dialog.addString("Nuclei IDs", "");
	Dialog.show();

	removed_nuclei = Dialog.getString(); 
	
	run("Replace/Remove Label(s)", "label(s)=" + removed_nuclei + " final=0");
	
	// Using the nuclei to get rid of any cell fragments
	run("8-bit");

	setThreshold(1, 255);
	run("Convert to Mask", "method=Li background=Dark black");
	rename("nuclei_mask");

	run("Morphological Reconstruction 3D", "marker=nuclei_mask mask=[microglia_binary] type=[By Dilation] connectivity=6");
	rename("microglia_singles");
	
	// Create and save a 3D-like image for microglia in red and nuclei in cyan
	run("Merge Channels...", "c1=microglia_singles c5=nuclei_mask create keep ignore");
	run("Z Project...", "projection=[Standard Deviation]");
	run("RGB Color");
	saveAs("Tiff", Output_folder + "/Max_" + filename + "_3DinRed.tif");
	close();
	
	// Separating touching cells using 3D watershed and saving a 3D-like image with each cell having a distinct color
	selectWindow("nuclei_mask");
	run("3D Maxima Finder", "minimmum=0 radiusxy=1.50 radiusz=1.50 noise=100");
	rename("seeds");
	
	run("3D Watershed", "seeds_threshold=254 image_threshold=254 image=microglia_singles seeds=seeds radius=2");
	setVoxelSize(width, height, depth, unit);
	run("glasbey_on_dark");
	run("Region Boundaries Labeling");
	run("RGB Color");
	run("Z Project...", "projection=[Standard Deviation]");
	rename("microglia_RGB");
	
	selectWindow("nuclei_mask");
	setVoxelSize(width, height, depth, unit);
	run("Z Project...", "projection=[Standard Deviation]");
	run("blue_orange_icb");
	run("RGB Color");
	rename("nuclei_RGB");
	imageCalculator("Add create", "microglia_RGB", "nuclei_RGB");
	
	saveAs("Tiff", Output_folder + "/Max_" + filename + "_Segmented_Labels.tif");
	close();

	// Closing unnecessary images
	close("Composite");
	close("nuclei_mask");
	close("nuclei_RGB");
	close("watershed-bnd");
	close("seeds");
	close("Result-Opening");
	close("Result-Opening-lbl");
	close("Result of microglia_binary");
	close("DAPI_RAW");
	close("DAPI_binary");
	
	// Creating and saving a 3D-like image for microglia cells with 2D skeleton overlaid
	selectWindow("microglia_singles");
	setVoxelSize(width, height, depth, unit);
	run("Z Project...", "projection=[Max Intensity]");
	run("Skeletonize (2D/3D)");
	run("Dilate");
	run("RGB Color");
	rename("branches_RGB");
	imageCalculator("Add", "branches_RGB", "microglia_RGB");
	saveAs("Tiff", Output_folder + "/Max_" + filename + "_Segmented_Branches.tif");
	close();
	
	// Analyze 3D morphological features using MorphoLibJ
	selectWindow("watershed");
	setVoxelSize(width, height, depth, unit);
	run("Analyze Regions 3D", "voxel_count volume surface_area mean_breadth sphericity euler_number bounding_box centroid equivalent_ellipsoid ellipsoid_elongations max._inscribed surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
	number_of_labels = getValue("results.count");
	
	// Adding columns to save the skeletonization data
	selectWindow("watershed-morpho");
	Table.set("# Branches", 0, 0);
	Table.set("# Junctions", 0, 0); 
	Table.set("# End-point voxels", 0, 0); 
	Table.set("# Junction voxels", 0, 0); 
	Table.set("# Slab voxels", 0, 0); 
	Table.set("Average Branch Length", 0, 0); 
	Table.set("# Triple points", 0, 0); 
	Table.set("# Quadruple points", 0, 0); 
	Table.set("Maximum Branch Length", 0, 0); 
	Table.set("Longest Shortest Path", 0, 0);
	Table.set("spx", 0, 0); 
	Table.set("spy", 0, 0); 
	Table.set("spz", 0, 0);  
	
	
	headingsArray_skel = newArray("# Branches", "# Junctions", "# End-point voxels", "# Junction voxels", 
							"# Slab voxels", "Average Branch Length", "# Triple points", 
							"# Quadruple points", "Maximum Branch Length", "Longest Shortest Path", 
							"spx", "spy", "spz");
	
	// Isolating sigle cells to run the skeleton analysis independently for each cell
	// Here the single skeleton images are not saved. To save them, uncomment the commented commands.
	
	for(i = 1; i < number_of_labels + 1; i++) {
		selectWindow("watershed");
		run("Crop Label", "label=" + i + " border=1"); 
		setVoxelSize(width, height, depth, unit);
		//saveAs("Tiff", Output_folder + "/" + filename + "_Cell.ID" + i + "_Segmented.tif");
		
		run("Z Project...", "projection=[Max Intensity]");
		//saveAs("Tiff", Output_folder + "/Max_" + filename + "_Cell.ID" + i + "_Segmented.tif");
		//selectWindow(Cell_ID_Seg);
		run("Skeletonize (2D/3D)");
		Cell_ID_Seg = getTitle();
		//saveAs("Tiff", Output_folder + "/" + filename + "_Cell.ID" + i + "_Skeletonized.tif");
		
		run("Analyze Skeleton (2D/3D)", "prune=none calculate");
		
		selectWindow("Longest shortest paths");
		//saveAs("Tiff", Output_folder + "/Max_" + filename + "_Cell.ID" + i + "_LS path.tif");
		close();
		//selectWindow("Longest shortest paths");
		//close();
		
		selectWindow("Tagged skeleton");
		//run("Z Project...", "projection=[Max Intensity]");
		//saveAs("Tiff", Output_folder + "/Max_" + filename + "_Cell.ID" + i + "_Tagged skeleton.tif");
		//close();
		//selectWindow("Tagged skeleton");
		close();
		
		selectWindow(Cell_ID_Seg);
		close();
			
			// Copying results to the main results table
			for(j = 0; j < headingsArray_skel.length; j++) {
			selectWindow("Results");
			data = Table.get(headingsArray_skel[j], 0);
			selectWindow("watershed-morpho");
			Table.set(headingsArray_skel[j], i-1, data);
			Table.update;
			}				 			 	
	}
	
	// Saving the 3D z-stack of the fully segemented 3D image
	selectWindow("watershed");
	run("glasbey_on_dark");
	saveAs("Tiff", Output_folder + "/" + filename + "_Segmented_Labels_zStack.tif");
	
	// Saving the results table as csv file
	selectWindow("watershed-morpho");
	saveAs("Results", Output_folder + "/" + filename + "_Results.csv");
	close(filename + "_Results.csv");
	
	// Clearing results and close all images and tables
	run("Clear Results");
	close("Results");
	close("*");
	
	// Inform user about progress
	print("The analysis of image " + n+1 + " out of " + lengthOf(fileList) + " is done!");
}
