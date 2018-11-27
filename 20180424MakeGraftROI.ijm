//importing working directory
path=File.openAsRawString("/Users/ryutaro/Dropbox/path.txt");
print(path)

//setting working directory of Fiji
open(path + File.separator + "image.tif");
name = getTitle();
print(name);

//Making graft ROI
setTool("polygon");
waitForUser("Select graft area, then hit OK.");  
roiManager("Add")
roiManager("Save", path + File.separator + "ROI_graftarea.zip");
run("Close All");
open(path + File.separator + "image.tif");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Reset");
run("Select None");
saveAs("Tiff", path + File.separator + "graft.tif");

//print("quitting...");
//run("Quit");
print("quitting...");
eval("script", "System.exit(0);"); 