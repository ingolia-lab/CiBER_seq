// function to auto-segment nuclei and measure DAPI and GCN4

/*
 * Macro template to process multiple images in a folder
 */

#@ File (label = "~/Desktop/Zuriah_analyzed_that/Sum_z_project/", style = "directory") input
#@ File (label = "~/Desktop/Zuriah_analyzed_that/nuclei_dapi_gcn4/", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	run("Bio-Formats", "open=" + input + File.separator + file + " autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	//run("Bio-Formats", "open=" + input + File.separator + file);
	title = getTitle();
	run("Duplicate...", "duplicate");
	//figure out how to choose second channel for DAPI Stack.setPosition(channel, slice, frame) 
	Stack.setPosition(2, 1, 1)
	setAutoThreshold("Otsu dark");
	//run("Threshold...");
	run("Convert to Mask", "method=Otsu background=Dark calculate");
	//setThreshold(255, 255);
	run("Convert to Mask", "method=Otsu background=Dark calculate only black");
	run("Analyze Particles...", "size=0.5-Infinity show=Outlines display exclude clear include add slice");
	// xls_name=replace(title,".tif",".xls");
	temp = replace(file, ".tif", "-1.tif");
	selectWindow("Drawing of " + temp);
	close();
	close();
	close("Results");
	selectWindow(title);
	run("ROI Manager...");
	roiManager("Deselect");
	roiManager("Measure");
	//roiManager("Select", 0);
	//run("Select All");
	// select all points in manager
	saveAs("Results", output + File.separator + file + "nuclei.csv");
	//roiManager("Measure");
	//saveAs("Results", output + File.separator + file + "nuclei.csv");
	close("Results");
	selectWindow(title);
	//new code
	selectWindow(title);
	run("Split Channels");
	selectWindow("C1-" + title); 
	roiManager("Deselect"); 
	roiManager("Measure");
	//new code
	saveAs("Results", output + File.separator + file + "nuclei_GCN4.csv");
	roiManager("Deselect");
	roiManager("Delete");
	// Leave the print statements until things work, then remove them.
	print("Processing: " + input + File.separator + file);
	print("Saving to: " + output);
}

// See also Process_Folder.py for a version of this code
// in the Python scripting language.
processFolder(input);

