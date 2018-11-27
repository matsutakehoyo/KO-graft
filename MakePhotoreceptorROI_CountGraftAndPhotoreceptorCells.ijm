roiManager("Reset");
run("Clear Results");

//---------!!You need to name the parent folder with desired condition value!!---------------
//importing working directory
path=File.openAsRawString("/Users/ryutaro/Dropbox/path.txt");
print(path)

setBatchMode(true);

//Process DAPI image
open(path + File.separator + "graft.tif");
name = getTitle();
print(name);
run("Duplicate...", "duplicate channels=1");
saveAs("Tiff", path + File.separator +  "dapi_graft.tif");
dapiZ=getImageID();
run("Duplicate...", "duplicate");
//run("Subtact background...")
run("Subtract Background...", "rolling=18");
//run("Bandpass filter...)
run("Bandpass Filter...", "filter_large=60 filter_small=6 suppress=None tolerance=5 process");
dapiZbandpass=getImageID();

//create empty image
selectImage(dapiZbandpass);
run("Duplicate...", "duplicate");
run("Select All");
run("Clear", "slice");
dapiEmpty=getImageID();
wait(800);

//Find maxima for segmentation
selectImage(dapiZbandpass);
run("Find Maxima...", "noise=3 output=[Segmented Particles]");
run("Analyze Particles...", "size=0-infinity clear add");
//roiManager("Save", path + File.separator + "PreMaxROI.zip");
nROI = roiManager("count");
print(nROI);

iniT =8;

//Loop for Every ROI and Threshold Skeletonize
for(i=0; i<nROI; i++){
	print(i);
	selectImage(dapiZbandpass);
	run("Duplicate...", "duplicate");
	temp=getImageID();
	roiManager("select", i);
	run("Clear Outside");
	getStatistics(area, mean, min, max, std, histogram);
	tol_i=(max - min)/1.5;
	if(tol_i<iniT){
		tol_i=iniT;
	}
	print(tol_i);
	run("Select None");
	run("Find Maxima...", "noise=&tol_i output=[Maxima Within Tolerance]");
	temp2=getImageID();
	run("Select All");
	getStatistics(area, mean);
	//print(mean);
	run("Select None");
	if(mean>254){
	run("Invert");
	}
	//run("Auto Threshold", "method=Otsu white");
	//run("Options...", "iterations=1 count=1 black do=Skeletonize");
	//run("Options...", "iterations=1 count=1 black do=Dilate");
	imageCalculator("Add add", dapiEmpty, temp2);
	selectImage(temp);
	run("Close");
	selectImage(temp2);
	run("Close");
	//Clear memory
	//call("java.lang.System.gc");
}
selectImage(dapiEmpty);
run("Adjustable Watershed", "tolerance=0.5");
saveAs("Tiff", path + File.separator +  "dapi_graft_watershed.tif");

//Make DAPI roi
run("Analyze Particles...", "clear add slice");
roiManager("Save", path + File.separator + "ROIset_dapi_graft.zip");
run("Close");
//Measure
setBatchMode("exit & display");
open(path + File.separator + "dapi_graft.tif");
roiManager("Show None");
roiManager("Show All");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction limit display scientific add nan redirect=None decimal=9");
roiManager("Measure");
saveAs("Results", path + File.separator +  "DAPI_graft.csv");
run("Clear Results");
roiManager("Reset");


//Making Photoreceptor ROI
open(path + File.separator + "graft.tif");
run("Duplicate...", "duplicate channels=2");
GFP=getImageID();
run("Smooth");
run("Smooth");
run("Smooth");
run("Smooth");
run("Smooth");
run("Smooth");
run("Smooth");
run("Smooth");
run("Smooth");
run("Smooth");
run("Auto Threshold", "method=Otsu white");
run("Options...", "iterations=2 count=1 black do=Dilate");
run("Analyze Particles...", "clear add slice");
roiManager("Combine");
roiManager("Delete");
roiManager("Add");
roiManager("Show None");
roiManager("Show All");
roiManager("Save", path + File.separator + "ROI_photoreceptorarea.zip");
//Make dapi photoreceptor ROI
open(path + File.separator + "dapi_graft_watershed.tif");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Reset");
run("Select None");
saveAs("Tiff", path + File.separator + "photoreceptor_watershed.tif");
run("Analyze Particles...", "clear add slice");
roiManager("Save", path + File.separator + "ROI_photoreceptordapi.zip");
//Measure
setBatchMode("exit & display");
open(path + File.separator + "dapi_graft.tif");
roiManager("Show None");
roiManager("Show All");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction limit display scientific add nan redirect=None decimal=9");
roiManager("Measure");
saveAs("Results", path + File.separator +  "DAPI_photoreceptor.csv");
run("Clear Results");
roiManager("Reset");


//Clear memory
call("java.lang.System.gc");

//print("quitting...");
//run("Quit");
print("quitting...");
eval("script", "System.exit(0);"); 
