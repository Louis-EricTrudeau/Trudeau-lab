Dialog.createNonBlocking("Count microglia in brain slices");
 
Dialog.addDirectory("Select your directory ", "/");

Dialog.addCheckbox("See process?", false);

Dialog.show();

masterDir = Dialog.getString();
seeProcess = Dialog.getCheckbox();

close("*");

if (seeProcess == 1) {} else {
	setBatchMode("hide");
}

dirList = getFileList(masterDir);
for (f = 0; f < dirList.length; f++){
	parentFile = masterDir + dirList[f];
	if(lengthOf(dirList[f]) < 10) { 

		dataDir = masterDir + dirList[f] + "/data/";
		File.makeDirectory(dataDir);
		
		Iba1SomaDir = dataDir + "Iba1 soma/";
		File.makeDirectory(Iba1SomaDir);
		
		input_images = masterDir + dirList[f] + "/Images/"; 
		list = getFileList(input_images);
		
		for (i = 0; i < list.length; i++){	
			if(lengthOf(list[i]) < 12) {
				file = input_images + list[i];
				
				run("Bio-Formats Importer", "open=[" + file + "] color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1 z_begin=1 z_end=11 z_step=1");
				Iba1 = getTitle(); 
				
				run("Z Project...", "projection=[Average Intensity]"); close(Iba1); projectSoma = getTitle(); 
				resetMinAndMax(); 
				run("Set Scale...", "distance=4.8309 known=1 unit=micron");
				run("Pseudo flat field correction", "blurring=50 hide");
				
				run("Subtract Background...", "rolling=25");
				run("Enhance Contrast...", "saturated=0.35");
				run("Apply LUT");
				setThreshold(21000, 65535);  
				run("Convert to Mask");
				run("Fill Holes");
				run("Erode");
				run("Open");
				run("Analyze Particles...", "size=10-Infinity show=[Masks]");
				run("Gaussian Blur...", "sigma=2");
				run("Convert to Mask");
				run("Open");
				run("Dilate");
				run("Analyze Particles...", "size=10-Infinity show=[Overlay Masks] add"); 
				close("*");
				
				run("Bio-Formats Importer", "open=[" + file + "] color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1 z_begin=1 z_end=11 z_step=1");
				run("Z Project...", "projection=[Max Intensity]");
				resetMinAndMax();
				run("Set Scale...", "distance=4.8309 known=1 unit=micron");
				
				roiManager("Deselect");
				run("Set Measurements...", "area mean perimeter shape redirect=None decimal=3");
				roiManager("multi-measure measure_all");
				
				roiManager("Deselect"); roiManager("Delete");
				close("*");
				
				selectWindow("Results");
				saveAs("Results", Iba1SomaDir + list[i] + ".csv");
				close("Results");
			}
			else { 
				i++;
			}
		}
	}
}