Dialog.createNonBlocking("Measure TH and DAT signal in brain slices");
 
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
	if(lengthOf(dirList[f]) < 8) { 
		input_images = masterDir + dirList[f] + "/Images/"; 
		list = getFileList(input_images);
		
		function newImage() {
				resetMinAndMax();
				run("Set Measurements...", "mean area_fraction limit redirect=None decimal=3");
		}
		
		for (i = 0; i < list.length; i++){	
			if(lengthOf(list[i]) < 12) {
				file = input_images + list[i];
				// Measure TH intensity
				run("Bio-Formats Importer", "open=["+file+"]"); TH = getTitle(); title = split(TH, "/"); title = replace(title[2], ".oif", "");
				newImage();
				run("Enhance Contrast...", "saturated=0.5"); run("Apply LUT", "slice");
				run("Convoluted Background Subtraction", "convolution=Gaussian radius=20 slice");
				setThreshold(18500, 65535);
				run("Create Selection");
				roiManager("Add");
				close("*");
				run("Bio-Formats Importer", "open=["+file+"]"); TH = getTitle(); title = split(TH, "/"); title = replace(title[2], ".oif", "");
				newImage();
				roiManager("Select", 0);
				run("Make Inverse");
				run("Subtract...", "value=65535 slice");
				roiManager("Deselect"); roiManager("Delete"); run("Select None");
				setThreshold(10, 65535); 
				run("Measure");
				setResult("Slice", nResults - 1, title); 
				IJ.renameResults("Summary");
				
				// Measure tdTomato intensity
				run("Next Slice [>]"); tdTomato = getTitle(); 
				newImage();
				run("Enhance Contrast...", "saturated=0.5"); run("Apply LUT", "slice");
				run("Convoluted Background Subtraction", "convolution=Gaussian radius=20 slice");
				setThreshold(18500, 65535);
				run("Create Selection");
				roiManager("Add");
				close("*");
				run("Bio-Formats Importer", "open=["+file+"]"); tdTomato = getTitle(); title = split(tdTomato, "/"); title = replace(title[2], ".oif", "");
				run("Next Slice [>]");
				newImage();
				roiManager("Select", 0);
				run("Make Inverse");
				run("Subtract...", "value=65535 slice");
				roiManager("Deselect"); roiManager("Delete"); run("Select None");
				setThreshold(10, 65535); 
				run("Measure"); tdTomatoIntensity = getResult("Mean", 0); tdTomatoArea = getResult("%Area", 0);
				close("Results");
				selectWindow("Summary"); IJ.renameResults("Results");
				setResult("tdTomato", nResults - 1, tdTomatoIntensity); 
				setResult("tdTomato area", nResults - 1, tdTomatoArea); 
				IJ.renameResults("Summary");
				
				// Measure DAT intensity
				run("Next Slice [>]"); DAT = getTitle(); 
				newImage();
				run("Enhance Contrast...", "saturated=0.5"); run("Apply LUT", "slice");
				run("Convoluted Background Subtraction", "convolution=Gaussian radius=20 slice");
				setThreshold(17000, 65535);
				run("Create Selection");
				roiManager("Add");
				close("*");
				run("Bio-Formats Importer", "open=["+file+"]"); DAT = getTitle(); title = split(DAT, "/"); title = replace(title[2], ".oif", "");
				run("Next Slice [>]");
				newImage();
				roiManager("Select", 0);
				run("Make Inverse");
				run("Subtract...", "value=65535 slice");
				roiManager("Deselect"); roiManager("Delete"); run("Select None");
				setThreshold(10, 65535); 
				run("Measure"); DATIntensity = getResult("Mean", 0); DATArea = getResult("%Area", 0);
				close("Results");
				selectWindow("Summary"); IJ.renameResults("Results");
				setResult("DAT", nResults - 1, DATIntensity); 
				setResult("DAT area", nResults - 1, DATArea); 
		
				close("*");
			} else { 
				i++;
			}
		}
		
		saveAs("Results", parentFile + "intensity.csv");
		close("Results");
	}
}