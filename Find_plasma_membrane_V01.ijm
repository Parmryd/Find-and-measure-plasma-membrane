/* Cell membrane - identify plasma membrane (PM) (or similar) ROI for two channels.  
Problems?       Contact Ingela Parmryd <ingela.parmryd@gu.se>
Macro used in https://doi.org/10.1016/j.bbamem.2022.184094

The macro expects two images(channels) - could be a Z-series of equal size.
Once one channel is loaded the macro tries to open the 2nd, will succeed if the names are similar and
a differences between the paired names is built in, i.g. 
  "Ch_longer"   &   "Ch_shorter" (relevant for GP-analysis)
  "483"    &    "442"
  "ch00"   & "ch01"
  NameLongerRef   NameShorterRef 
  Other name pairs can be used - change lines 71-72.

The images are aligned using Translation aligment that is built on StackReg 
These two plugins need to be installed for the macro to work and
can be found at https://github.com/fiji-BIG/StackReg and http://punias.free.fr/ImageJ/Translation_Alignment.java
Please cite
P. Thévenaz, U.E. Ruttimann, M. Unser, "A Pyramid Approach to Subpixel Registration Based on Intensity," 
IEEE Transactions on Image Processing, vol. 7, no. 1, pp. 27-41, January 1998. 
Other relevant publications on alignment are available at http://bigwww.epfl.ch/publications/

Find_plasma membrane - requires roundish cells/objects.
Initially a few sequential points are marked, low precision is sufficient.
The ROI is automatically generated
A search around each initial point for the most intense pixel is performed using one or three options:
(0) in a circular area round the initial point
(1) along a line perpendicular the plasma membrane
(2) along line from the center of the cell (the default setting)

Additional points are interpolated with selectable spacing 
and the local search for intensity maxima is repeated for the interpolated points.

ROI names in ROI Manager & their position(index) in ROI Manager when all ROIs exist
Crop          			CropROI - area extracted from original image 
Initial					InitialROI    - a few manually entered points
Centre					CentreROI  - the cell centre
Initial Mx				InitialMxROI - automatically optimised initial points
Initial Mx Edt  	 Initial Max Edt  - optional manual edit
Mmb Intp				MmbIntROI    - interpolated points
Mmb Intp Mx				MmbIntpMxROI - optimised interpolated points
Mmb Intp Mx Edt  	 Mmb Intp Mx Edt - optional manual edit

When more than one cell is cut from the same stack or image set, 
their origin is marked and retained using an overlay. 
Th cropped images are saved and when repeated measurements are performed from the same orignal images,
the name of the saved data is annotated - cell1, cell2 etc. 

cut down version of Generate_plasma_membrane_ROI_05f.ijm
*/

macroversion="Plasma_membrane_ROI_V01.ijm";// name of version, print with results
// 04c  turn off linear interpolation
// 04e optimise for for single imges
// 05  remove backup alignment, add 13ccutdown01e

run("Close All");
print("\\Clear");
run("Point Tool...", "type=Dot color=Orange size=Tiny counter=0");// point tool settings
// Search for maximal intensity.
checkInit=2; // Initial points: 0-circular area, 1-perpendicular, 2-from centre.
checkIntrp=2;// Interpolated points: 0-circular area, 1-perpendicular, 2-from centre.
  rngInit=3;  // Initial points search distance/radius, in pixels.
  rngIntrp=3; // Interploated points search distance/radius, in pixels.
SliceSearchUse=1; // Which image to use for PM search: 1-channel 1, 2-channel 2, 3-combined image.
sig=0.8; // Sigma for Gaussian blur when smoothing image, in pixels.

spacing=4; // Spacing for interpolation, in pixels.
pxlSze=0.091;// Pixel size in micrometers.
defaultSlc=0; // Default slice, 0 will select the middle slice from stack.
useSlc=defaultSlc;
NameLongerRef="LLong";// Include in name of Longer wavelength - you can alter the text.
NameShorterRef="LShrt";// Include in name of Shorter wavelength - you can alter the text.

run("Conversions...", " ");// no change to values when altering image type
roiManager("reset");// remove any existing ROIs
print("\\Clear");// empty Log window
print("macro ",macroversion);
if (nImages>0) run("Close All");// close all images
// END OF Default & settings

// Open Ch_longer image - longer wavelength
waitForUser("OPEN   Ch_longer image","....then 'OK'");// choose image
run("Select None");

// Find second channel   based on name of first
ImageFrom=getDirectory("image"); // use to look for second image and for saving
getRawStatistics(nPixels, mean, minA, maxA, std, histogram);
	setMinAndMax(minA, maxA/2);// contrast stretch display
OrigRtID=getImageID;
OrigName=getTitle(); // name of first image
// alter name to automatically load second channel  
OrigName1=replace(OrigName, "Ch_longer", "Ch_shorter");// now look for "Ch_shorter" use to open 2nd channel
	//print("look for ",OrigName1);
OrigName1=replace(OrigName1, "483", "442");// 2nd change to name
	//print("look for ",OrigName1);
OrigName1=replace(OrigName1, "ch00", "ch01");// 3rd change to name
	//print("look for ",OrigName1);
OrigName1=replace(OrigName1, NameLongerRef,NameShorterRef);// 4th change to name
	//print("original ",OrigName);// 
	//print("look for ",OrigName1);
	//print("name alter done");// name of second image
ImageFrom=File.directory;
	print("imagefrom ",ImageFrom);
rename("Right");
getDimensions(width, height, channels, slices, frames);
defaultSlc=minOf(1,round(slices/2));// middle slice  
print("defaultSlc",defaultSlc);
setSlice(defaultSlc);
useSlc=defaultSlc;
setSlice(defaultSlc);
if (slices>1 && defaultSlc==0) setSlice(round(slices/2));// use middle slice if no default
if (defaultSlc!=0 && defaultSlc<=slices) setSlice(defaultSlc); // use default
	print("slices",slices);
	if (slices>1) waitForUser("select slice: then OK"); // SELECT/APPROVE SLICE
	useSlc=getSliceNumber();// whatever is the current
//if (slices>1) useSlc=getSliceNumber();// whatever is the current
	//print("line 129 useSlc",useSlc);
if (slices>=1) run("Duplicate...", " ");// copy a single slice 
	run("Remove Overlay");
//print("slice for right image ",useSlc);
rename("Right single");
run("Select None");// remove any ROI

// lEFT image - its name expected to be derivable from first channel's (lines 68-73)
open(ImageFrom+OrigName1);
	if(nImages>3) { // second channel was not automatically found
		waitForUser("Problem: 2nd image missing", "Open manually......OK");
	}
rename("Left");
OrigLfID=getImageID;
getDimensions(width, height, channels, slices, frames);
setSlice(useSlc); // same slice as Right stack
	print("useSlc",useSlc);
	useSlc=getSliceNumber();
	print("slice used for left image",useSlc);
	// save chosen slice in an a selection
	SliceArray=newArray(1);
	SliceArray[0]=useSlc;
	print(SliceArray.length);
		makeSelection("point",SliceArray,SliceArray);
		setSelectionName("Slice"+useSlc);// crop original images
		roiManager("add");
		run("Select None");// no selection
//if (slices>1) setSlice(useSlc); // same as Right
run("Duplicate...", " ");// single slice 
run("Select None");
rename("Left single");

// MAKE Stack
run("Images to Stack", "name=Stack title=single use");//stack with 2 images (right single & left single)
setSlice(1);// right image
run("Select None");
getRawStatistics(nPixels, mean, minDispl, maxDispl);
	//print("read",minDispl, maxDispl);
maxDispl=maxDispl-minDispl-0.6*(maxDispl-minDispl);// stretch display
setMinAndMax(minDispl, maxDispl);
	//print("display",minDispl, maxDispl);
//waitForUser("select cell - ROI or none:  then OK");// select area - avoid areas that could mess up registration
setTool("rectangle");
// select area - try to avoid areas that might mess up registration 
waitForUser("       Choose single Cell"," (i) select area ............. then OK  \n or \n (ii) to use whole image ......just OK");
//if (selectionType>=0) run("Duplicate...", "duplicate");// crop stack if ROI;
CropROI=-1; // negative value flags absence of Cropped Area
if (selectionType>=0) {	// if there is a selection	
	run("Add Selection...");// add selection as an overlay - finally saved to show cells already used				
	setSelectionName("CropROI");// crop original images
	//roiManager("reset");// remove any existing ROIs
	roiManager("add");
	CropROI=roiManager("count")-1;
	run("Crop");  // crop stack of two images
} // if selection

// REGISTRATION 
run("StackReg ", "transformation=Translation");// align - rotation and translation
run("Translation Alignment"); // cut down version of StackReg plugin, additionally reports offset
//	xmove= Ext.getShiftValueX;
//	ymove= Ext.getShiftValueY;
//print("Registration  xmove",xmove,"ymove",ymove);
//makePoint(xmove,ymove); // make point selection that holds the alignment correction
//	setSelectionName("Align"); 
//	roiManager("Add");
//		run("Select None");
rename("Stack");// aligned images
setSlice(2);
run("Add Slice");// expand stack to hold a 3rd (combined image) available for membrane tracing
// Combine both images, an option for search for membrane
for (j=1;j<=2;j++) {  // adjust intensity - make similar-equal weight combined image
	selectImage("Stack");
	setSlice(j);
	run("Duplicate...", "use");
	if(j==1) rename("Right_Algn");
	if(j==2) rename("Left_Algn");
	run("Duplicate...", " ");
	// equalize intensities 
	if(j==1) rename("Right_AlgnScl");
	if(j==2) rename("Left_AlgnScl");
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	scl=30000/max;// max set to 30000
	run("Multiply...", "value=scl");
}
// combine the two channels
imageCalculator("Add create 32-bit", "Right_AlgnScl","Left_AlgnScl");// sum images
//imageCalculator("Add create 32-bit", "Right single","Left single");
run("Gaussian Blur...", "sigma=sig"); // smoothing
rename("CombAlignScl_Gau");
run("Select None");
run("Copy");
selectImage("Stack");
setSlice(3);
selectImage("CombAlignScl_Gau");
getLocationAndSize(xWin, yWin, widthWin, heightWin);
setLocation(10, 10, widthWin*2, heightWin*2);// zoom image and place in top left of screen

// Images combined        CHECK ALIGNMENT
run("Merge Channels...", "c1=Right_AlgnScl c2=Left_AlgnScl create");// make composite image
// merge also deletes the images - odd
// alter B & C
for (p=1;p<=2;p++) {
	setSlice(p);
	getRawStatistics(nPixels, meantmp, mintmp, maxtmp, std, histogram);
	setMinAndMax(mintmp, maxtmp);
}
Ht=getHeight();// size of image
Wd=getWidth();
setLocation(20, 20,Wd*2+54 ,Ht*2+82);// Zoom, note additions- for area of window incl slide bar
run("Cascade");// 
alignOK=getBoolean("Are Images correctly Aligned");// user prompt
//selectImage("Right_Algn");// right image - better image
if (alignOK==0) { // NOT ALIGNED - 
	print("not aligned   will quit");
	wait(1000);
	exit("Not aligned");
} // align
close();// delete composite image
	selectImage("CombAlignScl_Gau");
	print("Auto Align OK");
	run("Copy");
	selectImage("Stack");
	setSlice(3);
	run("Paste");
	run("Select None");
	Ht=getHeight();// size of image
	Wd=getWidth();
	setLocation(20, 20,Wd*2+54 ,Ht*2+82);// Zoom, note additions- for area of window incl slide bar	
	print("madeit");
	//selectImage("CombAlignScl_Gau");
	//close();
selectImage("Stack");	
setSlice(SliceSearchUse);// slice to search
run("Duplicate...", "use");
rename("SearchThis");
	if(sig>0) run("Gaussian Blur...", "sigma=sig");
Httemp=getHeight();// size of image
Wdtemp=getWidth();
//print("image",Wdtemp,Httemp);
setLocation(20, 20,Wdtemp*2+54 ,Httemp*2+60);// note additions-added area of window ZOOM 2
run("Select None");
getRawStatistics(nPixels, mean, min, max, std, histogram);
setMinAndMax(min, max*1.1);// alter display
setTool("polyline");
	getRawStatistics(nPixelsComb, meanComb, minComb, maxComb);
	maxComb=(maxComb-minComb)*0.90 +  minComb;
		//print(minComb,maxComb);
	setMinAndMax(minComb,maxComb);

// MARK INITIAL POINTs   precision not critical
selectImage("SearchThis");
run("Point Tool...", "type=Dot color=Cyan size=Tiny");
waitForUser("Membrane:  mark Initial points.....  then 'OK'");// FIRST ROI by hand...............
setSelectionName("Initial");
	// round initial coordinates to integer values
	getSelectionCoordinates(xInitArray, yInitArray);
	for(r=0;r<xInitArray.length;r++){
		xInitArray[r]=round(xInitArray[r]);	
		yInitArray[r]=round(yInitArray[r]);	
	}
	//Array.show(xInitArray);
roiManager("add");
InitialROI=roiManager("count")-1;

//   find Centre find Centre find Centre find Centre find Centre find Centre 
run("Fit Circle"); 
getBoundingRect(x, y, width, height);
cenX=x+width/2;
cenY=y+height/2;
	drawOval(cenX-1, cenY-1, 2,2);
	centXarray=newArray(1); 
		centXarray[0]=cenX;
	centYarray=newArray(1); 
		centYarray[0]=cenY;
	makeSelection(10, centXarray, centYarray);// point selection, central location
	setSelectionName("Centre");
	roiManager("add");
		CentreROI=roiManager("count")-1;

//     Optimise initial points
roiManager("select",InitialROI);
getSelectionCoordinates(xpt, ypt);// make arrays with x and y coordinate of marked points
// Use the  combined aligned smoothed image to optimise manually marked positions
// 3 search options    circle, perpendicular line, line from centre
rng=rngInit;// range for initial point search
roiManager("select", InitialROI);
getSelectionCoordinates(xpt, ypt);
check=checkInit;
//print("intial check");
run("Remove Overlay");
for (t=0;t<=xpt.length-1;t=t+1) { 
	//print("initial",t,"  check ",check);
	xat=xpt[t];// dot x pos
	yat=ypt[t];// dot y pos
	// find max I in circular area around initial point
	vOrig=getPixel(xat,yat);// read intensity value
	vMax=vOrig;// look for max intensity - first estimate
	//print (t,vOrig);
	//rng=5;// distance from marked points in pixels
	
	if (check==0) { // circular area around point
	  for (x=xat-rng;x<=xat+rng;x=x+1){ // vary x
		//print("vmax ",vMax);
		for (y=yat-rng;y<=yat+rng;y=y+1){ // vary y
			v=getPixel(x,y);
			//print(t,x,y,vmax);
			if (v>vMax) { // if a more intense value found
				vMax=v;
				xmax=x;ymax=y;// update max
				xpt[t]=x;ypt[t]=y; // save most intense pos			
			}// if v>max
			//print(t,x,y,vOrig,vMax);
		}// y	vary y pos
	  }// x vary pos
	} // if check==0 circular area
	
	if (check==1) { // perpendicular line
		if (t<xpt.length-1) {
			xat1=xpt[t+1];// next dot
			yat1=ypt[t+1];
		} // if t<xpt.elngth-1
		if (t==xpt.length-1) { // last point in list
			xat1=xpt[t-1]; yat1=ypt[t-1];// previous dot
		} // last point
		//print("Perpendicular ",xat,yat," to ",xat1,yat1);
		// find gradient of midline
		xgrad=xat-xat1; ygrad=yat-yat1;
		//print("gradients ",xgrad,ygrad);
		// single pixel steps
		scl=xgrad;
			if (abs(xgrad)<abs(ygrad)) scl=ygrad; //
		//scl=maxOf(abs(xgrad),abs(ygrad)); 
		//print(" scl ",scl);
			xgradUse=xgrad/abs(scl);
			ygradUse=ygrad/abs(scl);
			//print("grad scl ",xgradUse,ygradUse);
			//Overlay.drawLine(xat, yat, xat+ygradUse*30,yat+xgradUse*30);
			run("Overlay Options...", "stroke=red width=0 set");
			//Overlay.drawLine(xat-ygradUse*rng,yat+xgradUse*rng, xat+ygradUse*rng,yat-xgradUse*rng);
			//Overlay.show;
			// check for most intense pixel
			xst=xat-ygradUse*rng; // first location
			yst=yat+xgradUse*rng;
			mxVal=0;// max val
			run("Overlay Options...", "stroke=white width=0 set");
			for (l=0;l<rng*2;l=l+1) {
				x=xst+ygradUse*l;
				y=yst-xgradUse*l;
				v=getPixel(x,y);
					//Overlay.drawRect(x, y, 0.5, 0.5);
			//print(t,x,y,vmax);
		   	 if (v>vMax) {
				vMax=v;
				xmax=x; ymax=y;
				xpt[t]=x; ypt[t]=y;	// save most intense		
			 }// if
			} // l loop	
	} // check ==1   perpendicular check line

	if (check==2) { // along line to centre   cenX,cenY, 
		//print(t," From Centre ",xat,yat," to ",cenX,cenY);
		run("Overlay Options...", "stroke=red width=0 set");
		Overlay.drawLine(xat,yat,cenX,cenY);
		Overlay.show;
		// find gradient of midline
		xgrad=cenX-xat;
		ygrad=yat-cenY;
		//print("gradients ",xgrad,ygrad);
		// single pixel steps
		scl=xgrad;if (abs(xgrad)<abs(ygrad)) scl=ygrad; //max of X or Y disregarding sign
		//print(" scl ",scl);
			//xgradUse=xgrad/abs(scl);// perpendicular
			//ygradUse=ygrad/abs(scl);
			xgradUse=xgrad/abs(scl);
			ygradUse=ygrad/abs(scl);
				//print("rescale gradient ",scl,"xgraduse",xgradUse,"ygraduse",ygradUse);
			// check for most intense pixel
			xst=xat-xgradUse*rng; // first location
			yst=yat+ygradUse*rng;
			mxVal=0;// max val
			for (l=0;l<rng*2+1;l=l+1) {	// start outside cell
				x=round(xst+xgradUse*l);
				y=round(yst-ygradUse*l);
					//print("   test initial ",x,y);
					//run("Overlay Options...", "stroke=white width=0 set");
					//Overlay.drawRect(x, y, 0.5, 0.5);// OPTION TO DISPLAY
				v=getPixel(x,y);
				//print("   test initial ",x,y,"  v",v);
				run("Specify...", "width=1 height=1 x=x y=y slice=1");
					run("Add Selection...");
			//print(t,x,y,vmax);
			  if (v>vMax) {
				vMax=v;
				xmax=x; ymax=y;
				xpt[t]=x; ypt[t]=y;	// save most intense		
			  } // if
			} // l loop	
			//Overlay.show;   // show test points
	} // check==2   along line to centre
	Overlay.show;   // show test points
	//waitForUser("check gradient");
	run("Select None");
} //t loop    OPTIMISE Initial MARKED POINTS  OPTIMISE   POINTS   OPTIMISE MARKED POINTS
makeSelection(6,xpt,ypt);// first entry is selection type 6- segmented line, 7 - freehand
roiManager("add");
run("Hide Overlay");
InitialMxROI=roiManager("count")-1;//
roiManager("select",InitialMxROI);
//setSlice(1);
roiManager("rename", "Initial Mx");
roiManager("select",InitialMxROI);
origN=xpt.length;// original number of points
EditYN=getBoolean("Manual Edit Points?");// Option to manually edit points
//Edit Points
//print("move yn ",moveYN);
Hv_InitialMaxEdt=0;// no edited Initial Mx ROI
if (EditYN==1) {
	print("edit intitial points");
	run("Hide Overlay");
	waitForUser("OPTION", "move/delete initial points.... then 'OK'");
		getSelectionCoordinates(xpntTMPArray, ypntTMPArray);
		makeSelection(6,xpntTMPArray, ypntTMPArray);
		roiManager("add");
	InitialMxEdtROI=roiManager("count")-1;//
	//setSlice(1);
	roiManager("select",InitialMxEdtROI);
	roiManager("rename", "Initial Max Edt");
	Hv_InitialMaxEdt=1; // edited version exists
} // end of edit membrane

// Interpolate Interpolate Interpolate Interpolate Interpolate Interpolate Interpolate
run("Fit Spline");
run("Interpolate", "interval=spacing");// spacing set at start of macro
waitForUser("check spline", "Spline fitted");
seltype=selectionType() ;
//print("spline seltype",seltype);
getSelectionCoordinates(xMaxInt, yMaxInt);
makeSelection(2, xMaxInt, yMaxInt);// polygon
setSelectionName("Mmb Intp");
	for(r=0;r<xMaxInt.length;r++){ // round coordinates
		xMaxInt[r]=round(xMaxInt[r]);	
		yMaxInt[r]=round(yMaxInt[r]);	
	}
	makeSelection(2, xMaxInt, yMaxInt);// polygon
    setSelectionName("Mmb Intp");
	roiManager("add");
	MmbIntROI=roiManager("count")-1;
	Array.show(xMaxInt);
	
// INTERPOLATION DONE   INTERPOLATION DONE INTERPOLATION DONE INTERPOLATION DONE INTERPOLATION DONE 
getSelectionCoordinates(xpt, ypt);
// optimise interpolated points
rng=rngIntrp;// select range for interpolated search
for (t=0;t<=xpt.length-1;t=t+1) { // for each point in interpolated ROI
	xat=xpt[t];yat=ypt[t];// dot's pos
	// find max I in circular area around each of the initial points
	vOrig=getPixel(xat,yat);
	vMax=vOrig;
	//print (t,vOrig);
	//rng=5;// distance from marked points in pixels	- for demo
	if (check==0) { // circular area
	//print(t,"search in circular area");
	for (x=xat-rng;x<=xat+rng;x=x+1){ // vary x
		//print("vmax ",vMax);
		for (y=yat-rng;y<=yat+rng;y=y+1){ // vary y
			v=getPixel(x,y);
			//print(t,x,y,vmax);
			if (v>vMax) {
				vMax=v;
				xmax=x; ymax=y;
				xpt[t]=x; ypt[t]=y;	// save most intense		
			}// if
			//print(t,x,y,vOrig,vMax);
		}// y	vary y area
	}// x vary x area
	} // if check==0 circular area
	
	if (check==1) { // perpendicular line
		print(t,"search perpendicular");
		if (t<xpt.length-1) {
			xat1=xpt[t+1];// next dot
			yat1=ypt[t+1];
		} // if t<xpt.elngth-1
		if (t==xpt.length-1) { // last point in list
			xat1=xpt[t-1];// next dot
			yat1=ypt[t-1];
		} // last point
		//print("Perpendicular ",xat,yat," to ",xat1,yat1);
		// find gradient of midline
		xgrad=xat-xat1;
		ygrad=yat-yat1;
		//print("gradients ",xgrad,ygrad);
		// single pixel steps
		scl=xgrad;if (abs(xgrad)<abs(ygrad)) scl=ygrad; //
		//scl=maxOf(abs(xgrad),abs(ygrad)); 
		//print(" scl ",scl);
			xgradUse=xgrad/abs(scl);
			ygradUse=ygrad/abs(scl);
			//print("grad scl ",xgradUse,ygradUse);
			//Overlay.drawLine(xat, yat, xat+ygradUse*30,yat+xgradUse*30);
			Overlay.drawLine(xat-ygradUse*rng,yat+xgradUse*rng, xat+ygradUse*rng,yat-xgradUse*rng);
			Overlay.show;
			// check for most intense pixel
			xst=xat-ygradUse*rng; // first location
			yst=yat+xgradUse*rng;
			mxVal=0;// max val
			for (l=0;l<rng*2;l=l+1) {
				x=xst+ygradUse*l;
				y=yst-xgradUse*l;
				v=getPixel(x,y);
			//print(t,x,y,vmax);
			if (v>vMax) {
				vMax=v;
				xmax=x; ymax=y;
				xpt[t]=x; ypt[t]=y;	// save most intense		
			}// if
		} // l loop	
	} // check =1   perpendicular check line
	// ADD FINAL 3x3 pixel check

	if (check==2) { // along line from centre   cenX,cenY, 
		//print(t," towards Centre ",xat,yat," to ",cenX,cenY);
		// find gradient of midline
		xgrad=cenX-xat;
		ygrad=yat-cenY;
		//print("gradients ",xgrad,ygrad);
		// single pixel steps
		scl=xgrad;if (abs(xgrad)<abs(ygrad)) scl=ygrad; //  larger of X or Y 
		//scl=maxOf(abs(xgrad),abs(ygrad)); 
		//print(" scl ",scl);
			xgradUse=xgrad/abs(scl);
			ygradUse=ygrad/abs(scl);
				//print("scale ",scl,"xgraduse",xgradUse,"ygraduse",ygradUse);
			// FIND for most intense pixel
			xst=round(xat-xgradUse*rng); // first location
			yst=round(yat+ygradUse*rng);
			mxVal=0;// max val
			for (l=0;l<rng*2;l=l+1) {	// start outside cell
				x=round(xst+xgradUse*l);
				y=round(yst-ygradUse*l);			
				v=getPixel(x,y);
				// show tested pixels
				//run("Specify...", "width=1 height=1 x=x y=y slice=1");
				//run("Add Selection...");//
			//print(t,x,y,vmax);
			  if (v>vMax) {
				vMax=v;
				xmax=x;ymax=y;
				xpt[t]=x; ypt[t]=y;	// save most intense		
			  } // if
			} // l loop	
	} // check ==2   along line to centre
} //t loop runs through points   
// HAVE INTERPOLATED MAX  HAVE INTERPOLATED MAX   HAVE INTERPOLATED MAX   HAVE INTERPOLATED MAX

makeSelection(6,xpt,ypt);// first entry is selection type 6- segmented line, 7 -freehand line, 3- freehand area
	roiManager("add");
MmbIntpMxROI=roiManager("count")-1;
	roiManager("select",MmbIntpMxROI);
roiManager("rename", "Mmb Intp Mx");
waitForUser("Max of Interpolated points done");
//make into points
		getSelectionCoordinates(xpntTMPArray, ypntTMPArray);
		//makeSelection(10,xpntTMPArray, ypntTMPArray);
		//run("Point Tool...", "type=Dot color=Red size=Tiny counter=0");
		
// Edit Points for PM  Edit Points for PM   Edit Points for PM     Edit Points for PM
run("Hide Overlay");
EditYN=getBoolean("Manual Edit Interp Points?");// Option to manually edit points
if (EditYN==1) {
	run("Hide Overlay");
	waitForUser("OPTION", "move/delete/add Interp points: then 'OK'");
		getSelectionCoordinates(xpntTMPArray, ypntTMPArray);
		makeSelection(6,xpntTMPArray, ypntTMPArray);
		roiManager("add");
		MmbIntpMxEdtROI=roiManager("count")-1;
	roiManager("select",MmbIntpMxEdtROI);
	roiManager("rename", "Mmb Intp Mx Edt");
} // end of edit membrane


// SAVE SAVE SAVE   images and ROIs

// save stack of aligned images
selectImage("Stack");
print(" ");
print("original name",OrigName);
//stackNamePart=replace(OrigName, "Ch_shorter.tif","Stack_Algnd-CellNN");// save stack as..., in same folder
//baseName=replace(OrigName, "Ch_shorter.tif","Stack_Algnd-CellNNN.tif");// save stack as..., in same folder
baseName=replace(OrigName, ".tif","Stack_Algnd-CellNNN.tif");// save stack as..., in same folder
print("baseName",baseName);

//replace(string, "NNN", 9);
//save(ImageFrom+stackName);
// check if file name "....cellNN.tif" already exists- more than one cell can be taken from the original image
print("  ");
ListofStacksArray=getFileList(ImageFrom);//
nameExists=0; //initial value
for (k=1;k<=9;k++) {
	posname=replace(baseName,"NNN",k);
	//print("checking posname",posname);
	nameExists=0;// posname already used
	for (h=0;h<ListofStacksArray.length;h++) {
		checkName=ListofStacksArray[h];
		//print(h,"    checkname",checkName);
		matchIndex=indexOf(checkName, posname);// does image exist
		if (matchIndex==0) nameExists=1; // name already used
		//print(k,h,"matchIndex",matchIndex);
	// find a CellNN variant that is not present
	} // h loop
	if ( nameExists==0) break;// all checked, the name is not in use 
} // k loop
print("use this name +1...  ",posname);
StackSaveAt=ImageFrom+posname;
	print("stack saved at ",StackSaveAt);   
	save(StackSaveAt);
// save ROIs
ROIsSaveAt=replace(StackSaveAt,"Stack","ROIs");
ROIsSaveAt=replace(ROIsSaveAt,".tif",".zip");
//ROIsSaveAt=ImageFrom+ROIsSaveAt;
	print("ROIs save at");
	print(ROIsSaveAt);
	roiManager("Deselect");// note, in the absence of a selection, the whole Manager is saved
	roiManager("save",ROIsSaveAt);

// Show "CropROI" on all images in stack as an overlay-show area previously analyzed
if(CropROI>=0) {
selectImage("Right");
roiManager("select", CropROI);// load first ROI
  for(j=1;j<=nSlices;j++) {
	setSlice(j);
	run("Add Selection...");// add selection as overlay
  }
run("Select None");
// save marked stack
StackSaveAt=ImageFrom+OrigName;
save(StackSaveAt);
}// if CropROI

selectWindow("Log");
print("Macro Finished");
