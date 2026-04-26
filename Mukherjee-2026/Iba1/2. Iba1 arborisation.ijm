Dialog.createNonBlocking("Detect microglia arborisation in brain slices");
 
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
		
		Iba1ArboDir = dataDir + "Iba1 arborisation/";
		File.makeDirectory(Iba1ArboDir);
		
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
				run("Pseudo flat field correction", "blurring=5 hide");
				
				run("Subtract Background...", "rolling=10");
				run("Enhance Contrast...", "saturated=0.35");
				run("Apply LUT");
				run("Despeckle");
				setThreshold(9000, 65535); 
				run("Analyze Particles...", "size=2-Infinity show=Masks");
				run("Invert LUT");
				run("Create Selection");
				roiManager("Add");
				close("*");
				
				run("Bio-Formats Importer", "open=[" + file + "] color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1 z_begin=1 z_end=11 z_step=1");
				run("Z Project...", "projection=[Max Intensity]");
				resetMinAndMax();
				run("Set Scale...", "distance=4.8309 known=1 unit=micron");
				
				roiManager("Select", 0);
				run("Set Measurements...", "area mean limit redirect=None decimal=3");
				run("Measure");
				setResult("Image", nResults - 1, list[i]);
				
				roiManager("Deselect"); roiManager("Delete");
				close("*");
			}
			else { 
				i++;
			}
		}
		selectWindow("Results");
		saveAs("Results", Iba1ArboDir + "Iba1 arborisation.csv");
		close("Results");
	}
}