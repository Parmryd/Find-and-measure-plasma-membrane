/* Measurements from an existing ROI and a pair of images
 Problems?       Contact Ingela Parmryd <ingela.parmryd@gu.se>
Macro used in https://doi.org/10.1016/j.bbamem.2022.184094
Measures intensity profiles from ROIs at a range of distances from membrane.
The range is settable, e.g.-3 to +3 pixels with minus distances inside the cell.
Input needed: (i) a stack with at least two images and (ii) a set of ROIs created by the related macro Find_plasma_membrane_V01.ijm
Measurements - values from two fluorescent channels are read, here used to calculate GP-values:

GP=(I(shorter wavlength range)-I(longer wavelength range))/(I (shorter wavelength range)+I(longer wavelength range) )

Overall measurements are made from the whole of each offset membrane ROI.
In addition each pixel in a single-pixel-wide ROI derived from the each ROI  is read sequentially and added to the Results.
ROIs can be created at a range of distances around the original ROI based on a Distance Map.
Version cut down from     TurnPmPointsIntoDistanceMap16d
  
  LUT  "Cent grey-8colramp" is needed - available along with this macro at GitHub
*/

//  Tidy Up
if (isOpen("Results")) { // close results window, if open
	selectWindow("Results");
	close("Results");
}
run("Close All");// any images
roiManager("reset");// close ROI Manager
print("\\Clear");
// settings
run("Colors...", "foreground=white background=black selection=red");
version="Measure_plasma_membraneV01.ijm";
roiManager("reset");
run("Options...", "iterations=1 count=1 black edm=32-bit");// DT 32bit
// symmetrical range of distances around initial PM - in pixels
uRnge=getNumber("range around PM (pixels)", 2); //
lRnge=0-uRnge;// inwards, same distance

/* 
Function, SelectNamedROI(SearchForROIname)
Selects a named ROI - if not found, ends without a selection.
Afterwards
use     roiManager("index")  to save the ROI's position in the ROI Manager
use     selectionType()  to check whether the ROI was found
*/
function SelectNamedROI(SearchForROIname) {
nRois=roiManager("count");// number in ROI Manager
 for (h=0;h<nRois;h++) {
	roiManager("select", h); // in order to read its name 
		if (selectionName == SearchForROIname) break
	run("Select None");	
 } // h loop  look for ROI
} // end of function SelectNamedROI()

waitForUser("Open Image...then.'OK'");
getDimensions(width, height, channels, slices, frames);
OrigWdth=getWidth();// size of image
OrigHt=getHeight();

// Automatically Load ROIs	
UseDirectory=getDirectory("image");
title=getTitle();
	ROItitle=replace(title, "Stack", "ROIs");// modify image name to corresponding ROIs
	ROItitle=replace(ROItitle, ".tif", ".zip"); // modify image name to corresponding ROIs
		openROIs=UseDirectory+ROItitle;
		open(openROIs);	// load ROIs
//if (roiManager("count")==0) waitForUser("Problem", "No ROIs Loaded");
if (roiManager("count")==0) waitForUser("Auto Load ROIs failed, load....then.'OK'");

// Background intensities for each channel - create a ROI manually
// Automatically set display to show background
run("Select None");// whole image for stats
getRawStatistics(nPixels, mean, min, maxAll, std, histArray);
run("Restore Selection");
maxArray=Array.findMaxima( histArray, 5);
maxN=maxArray[0];
//upperI=min+(maxN-min)*40; //10% of range
upperI=maxAll/8;// want to see background
	//print("display range",min,upperI);
setMinAndMax(min, upperI);
run("16 colors");// false colour scale
// Drawing TOOL- select one of the next 3 lines
	setTool("oval");
	//setTool("polygon"); 
	//setTool("rectangle");
waitForUser("  Background values", "Outline an Area.....then OK");
setSelectionName("Bkgrnd");
roiManager("Add");
setSlice(1);
	getRawStatistics(nPixels, BKgrnd442, min, max, std, histogram);
	BKgrnd442=round(BKgrnd442);
	run("Grays");
	setMinAndMax(0, maxAll);
setSlice(2);
	getRawStatistics(nPixelsBK, BKgrnd483, min, max, std, histogram);
	BKgrnd483=round(BKgrnd483);
	//print("bck",BKgrnd442,BKgrnd483);			
run("Select None");// remove background area

// find membrane ROI
// print("search for    Mmb Intp Mx",selectionType()); 
SelectNamedROI("Mmb Intp Mx Edt");// edited version, may not exist
	selType=selectionType(); // check on whether a ROI was found
if(selType<0) SelectNamedROI("Mmb Intp Mx"); //no edited ROI, look for undedited version
	selType=selectionType();	
if(selType==-1) waitForUser("Problem", "No MembraneROI found");// Error message
MembUseIndex=roiManager("index");// save location of ROI in Manager
nameOfROI=selectionName;
	//print("selectionName of membrane ROI",nameOfROI);
getSelectionCoordinates(xLineArray, yLineArray); // read coordinates into arrays

// Create Results Table
setResult("Notes", 0,"          Mean");
setResult("Notes", 1,"           Std");
setResult("Notes", 2,title);// name of image
setResult("Notes", 3,version);// name of macro
setResult("Notes", 4,"Background lower channel:");//background
setResult("Notes", 5,BKgrnd442);//background
setResult("Notes", 6,"Background upper channel:");//background
setResult("Notes", 7,BKgrnd483);//background
setResult("Notes", 8,"Background pixels: "+nPixelsBK);//background
	getVoxelSize(Pxlwidth, Pxlheight, Pxlepth, Pxlunit); // pixel size
setResult("Notes", 9,"Pixel size=0.091"); //Pxlunit +" " + Pxlheight);//pixel
setResult("Notes", 10,"ROI "+nameOfROI);
setResult("Notes", 11,"   User notes:");
setResult("Notes", 12,"**********");
setResult("Notes",13,"Date:YearMonthDay");
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
month=month+1;
date="  "+year+"_"+month+"_"+dayOfMonth;
setResult("Notes",14,date);
lastNotes=14; // location of last entry in "Note" column
 setResult("xposO", 0, "null");// dummy values - forces position of column
 setResult("yposO", 0, "null");
 setResult("xposO", 1, "null");
 setResult("yposO", 1, "null");
 for(l=lRnge; l<=uRnge; l=l+1) {  // GP column headings for Results
 		//print("GP column headings l",l);
 	GPlabel="GP_" + l;
	setResult(GPlabel, 0,20);
	setResult(GPlabel, 1,21);
 }

newImage("Untitled", "8-bit black", OrigWdth, OrigHt, 1);

//  FIND Centre ROI
SelectNamedROI("Centre");// find ROI named "Centre"
	CentreROIIndex=roiManager("index");	
getSelectionCoordinates(xpoints, ypoints);// just one point
	//print("centre",xpoints.length);
Xcent=xpoints[0]; Ycent=ypoints[0];// save location of centre

// fill gap in PM - fill with circle between end points
// draw lines from centre to each end
Xst=xLineArray[0]; Yst=yLineArray[0];
		Overlay.drawLine(Xst,Yst, Xcent, Ycent);
LastEntry=xLineArray.length-1;//
Xend=xLineArray[LastEntry]; Yend=yLineArray[LastEntry];
	Overlay.drawLine(Xend,Yend, Xcent, Ycent);
	Overlay.show;
//print("end",Xend,Yend,"start", Xst, Yst);
// Radii to each end point
radius1=sqrt(pow((Xcent-Xst),2)+pow((Ycent-Yst),2));// start pnt to centre
radius2=sqrt(pow((Xcent-Xend),2)+pow((Ycent-Yend),2));// end pnt to centre
// Gradient from centre to start and end points
XgradSt=Xst-Xcent;
YgradSt=Yst-Ycent;
XgradEnd=Xend-Xcent;
YgradEnd=Yend-Ycent;
	//print("start grad",XgradSt,YgradSt);
	//print("end grad",XgradEnd,YgradEnd);
run("Remove Overlay");
for(r=0;r<10;r++) {
	newGradX=XgradEnd-(XgradEnd-XgradSt)*r/10 ;  //update gradients
	newGradY=YgradEnd-(YgradEnd-YgradSt)*r/10 ; 
	//print(r,r/10,round(newGradX),round(newGradY));
		newX=Xcent+newGradX*1;
		newY=Ycent+newGradY*1;
		// find distance and rescale to radius
		l1=sqrt(pow((Xcent-newX),2)+pow((Ycent-newY),2));
		stretch=radius1/l1;
			//print("       L1",l1,"stretch",stretch);
		newX=Xcent+newGradX*stretch;
		newY=Ycent+newGradY*stretch;
		  Overlay.drawLine(newX,newY, Xcent, Ycent);
		  Overlay.show;
		 // wait(300);// slow down to see what happens
} //
// Fill gap in PN, new points ROI, at radius from centre
NnewPnts=9; // number to add, could be changed
  GapFillarrayX=newArray(NnewPnts+2);// arrays for new gap filling points
  GapFillarrayY=newArray(NnewPnts+2);
  DTzapArrayX=newArray(NnewPnts+3);// exclude unwanted part of DTmap
  DTzapArrayY=newArray(NnewPnts+3);//  need to expand GapFillROI & include centre
 // insert start and end points - from PM ROI, can later combine with
  GapFillarrayX[0]=Xend;// start of PM
  GapFillarrayY[0]=Yend;
  GapFillarrayX[NnewPnts+1]=Xst;//  end of PM
  GapFillarrayY[NnewPnts+1]=Yst;
radius1=sqrt(pow((Xcent-Xst),2)+pow((Ycent-Yst),2));// start pnt to centre
radius2=sqrt(pow((Xcent-Xend),2)+pow((Ycent-Yend),2));// end pnt to centre
//for (h=1;h<=NnewPnts;h=h+1) { // each point used to fill gap in PM
for (h=0;h<=NnewPnts+1;h=h+1) { // each point used to fill gap in PM	
	frac1=h/(NnewPnts+1);// weight given to start pnt
	frac2=(NnewPnts+1-h)/(NnewPnts+1); //weight given to end pnt
	//print("h",h,"   ",h/(NnewPnts+1),"  ",((NnewPnts+1-h)/(NnewPnts+1)));
	MidX=Xst*frac1+Xend*frac2;// pos on to first-last line
	MidY=Yst*frac1+Yend*frac2;//  pos on to first-last line
	//Overlay.drawEllipse(MidX-1, MidY-1, 2,2);// circle at each new point
	//Overlay.show;
	 LtoMid=sqrt(pow((Xcent-MidX),2)+pow((Ycent-MidY),2));// L new pnt to centre
	 radiusUse=radius1*frac1+radius2*frac2; //use weighted radius - departures from circularity
	    xgradX=(MidX-Xcent);
	    ygradY=(MidY-Ycent);
	MidXext=Xcent+(radiusUse/LtoMid)*xgradX;// extend line to weighted radius
	MidYext=Ycent+(radiusUse/LtoMid)*ygradY;
	ZapMidXext=Xcent+(1.5*radiusUse/LtoMid)*xgradX;// extend line to weighted radius and outwards
	ZapMidYext=Ycent+(1.5*radiusUse/LtoMid)*ygradY;// extend by 50%
		//print("LtoMid",LtoMid,"  ");	
		Overlay.drawEllipse(MidXext-1, MidYext-1, 2, 2);// mark point on gap PM
		Overlay.show;
	    	GapFillarrayX[h]=MidXext;// save in Gap Array
		    GapFillarrayY[h]=MidYext; 
		    DTzapArrayX[h]=ZapMidXext;// Zap DT array
		    DTzapArrayY[h]=ZapMidYext;
}// h loop - fill gap with a few points 

// make DTZap ROI - use later to void part of Distance Map arising from filled gap
// triangular ROI, based on centre and two end points
// 
DTzapArrayX[NnewPnts+2]=Xcent;// adding central xy  - array made 20 lines earlier        
DTzapArrayY[NnewPnts+2]=Ycent;
makeSelection(2,DTzapArrayX,DTzapArrayY);// type 2, a polygon
setSelectionName("DTzap");
	roiManager("Add");
DTzapIndex=roiManager("count")-1;//location in ROI Manager
// smooth Gap filling line segment
makeSelection(6, GapFillarrayX,GapFillarrayY);// new ROI from start,end & new pnts	
		//maybe make ROI   type 2 polygon instead of type 6			
run("Fit Spline");// smooth ROI
setSelectionName("GapSpline");
	roiManager("Add");
GapSplineIndex=roiManager("count")-1;//location in ROI Manager
IndexArray=newArray(MembUseIndex,GapSplineIndex);// PM and FillGap Indexes -an array
roiManager("select", IndexArray); // select two ROIs	
roiManager("Combine");// combine (OR) 2 arrays
setSelectionName("PM-GapFilled");
	roiManager("Add");
		PM_GapFilledIndex=roiManager("count")-1;
	setForegroundColor(255,255,255);
run("Fill", "slice");// with 255, makes binary outline
run("Select None");// turn off ROI to allow filling

setForegroundColor(255,255,255);
run("Fill Holes");//   Binary FILLED CELL   Binary FILLED CELL   Binary FILLED CELL 

run("Create Selection"); // from binary cell
	setSelectionName("Whole Cell");
	roiManager("add");
	WholeCellIndex=roiManager("count")-1;	
	
//                     Make Distance Map - from binary
newImage("Cell_Binary", "8-bit black", width,height, 1);// 
roiManager("Select", PM_GapFilledIndex); // PM as a line
run("Fill", "slice");
	run("Clear Outside");// binary with PM outline only
run("Select None");
run("Invert");
run("Options...", "iterations=1 count=1 black edm=32-bit");// 32 bit distance map
setForegroundColor(255, 255, 255);
setBackgroundColor(0, 0, 0);
run("Distance Map");// from PM inwards and outwards
roiManager("select",WholeCellIndex);
run("Multiply...", "value=-1");// makes sign of inner area negative
	run("Add Selection...");// overlay from outline ROI
	run("Overlay Options...", "stroke=white width=0 set apply");
run("Select None");
getRawStatistics(nPixels, mean, minDT, maxDT);// to set display range
setMinAndMax(-30, 30);// could alter based on previous line, or base on uRnge & lRnge (distances
// remove gap related part of DTmap
// want DTzap  ROI
//getRawStatistics(nPixels, mean, minDT, maxDT);
setForegroundColor(maxDT,maxDT,maxDT);// value used for fill - 
	//print("maxDT",maxDT,"maxFill",maxFill);
roiManager("select", DTzapIndex);
changeValues(minDT, maxDT,maxDT); // reset Zapped area
run("Select None");// remove selection, PM remains as outline
// PM as overlay in black
roiManager("select", DTzapIndex-1);// outline with gap
run("Overlay Options...", "stroke=black width=0 fill=black set apply");
//run("Add Selection...");// add selection as overlay
run("Select None"); 
//run("Cent grey-8colramp 08a");// JA LUT false colour & greyscale(-ve values) ADD TO ImageJ
run("16 colors");
setMinAndMax(-8, 8);
rename("DistanceMapPosInOut");	
run("Hide Overlay");
//                                     HAVE Zapped/edited DISTANCE MAP


// MEASUREMENTS for range of distances from outline - ROIs added to Manager
// Put profiles in Results window
// show ROIs at distances away from original ROI
selectImage("DistanceMapPosInOut");
getDimensions(width, height, channels, slices, frames);

print("\\Clear");
// rescale gradients
XgradSt=XgradSt/radius1;
YgradSt=YgradSt/radius1;
setBatchMode(1);
for(l=lRnge; l<=uRnge; l=l+1) { // Range of distances from original perimeter
//for(l=-3; l<=-2; l=l+1) { // distances from original ROI - perimeter
	//print("  ");
	//print("position in Distance Loop l",l);
	row=0;
	heading1="L442 "+l;
	columnNameSlice1=heading1;
		setResult(heading1, 0, "null");// blank values
		setResult(heading1, 1, "null");// blank values
	heading2="L483 "+l;	
	columnNameSlice2=heading2;
		setResult(heading2, 0, "null");// blank values
		setResult(heading2, 1, "null");// blank values
			// Results table present
	
	// Distance Map - select part by thresholding
	selectImage("DistanceMapPosInOut");	
	run("Select None");//
	run("Duplicate...", " ");
	rename("Distance" + l);
	DTuseImageID=getTitle();
	// make distance ROI
	setThreshold(l-0.5, l+0.495);// one pixel range
	     //setThreshold(-0.5, 0.5);
	//print("range",l-0.5,l+0.5);
		run("Convert to Mask");
	        run("Skeletonize"); // new in version 12L
		run("Create Selection");
		//run("Create Selection");// why a repeat line????
		roiManager("Add");
		Nrois=roiManager("count");
		DisRoiIndex=Nrois-1;
		roiManager("select",DisRoiIndex);
		roiManager("Rename", "Dis"+l);
		// work on binary image, make binary, check neighbors, find end(s), track from end
		run("Divide...", "value=255");
		run("Convolve...", "text1=[1 1 1\n1 1 1\n1 1 1\n]");// number of Neighbors for each pixel  
		setMinAndMax(0, 4); // display range
		//have image with pixel connection numbers
		// find start point in binary image = connections ==2
		// run along line from centre
	Xmean=Xcent; // does this help
	Ymean=Ycent;
		Xsearch2=Xmean;// centre
		Ysearch2=Ymean;
		fnd=0;// flag when end point found
		//for (s=lngth-8;s<lngth+8;s=s+0.5) {  // search for start point, value 2 pixel
		//print("    XgradSt",XgradSt,"YgradSt",YgradSt,"Xmean",Xmean,"Ymean",Ymean);
		
		  for (s=radius1-8;s<radius1+8;s=s+0.5) {  // search for start point, value 2 pixel
		  	//print("s",s);
			for (x=-1;x<=1;x++) {
				//print("s",s,"x",x);
				for (y=-1;y<=1;y++) {
					//print("s",s,"x",x,"y",y);
			//print("s",l,s);
			//xtry=round(x+Xmean+s*xgradX1);XgradSt
			//ytry=round(y+Ymean+s*ygradY1);
			xtry=round(x+Xmean+s*XgradSt);
			ytry=round(y+Ymean+s*YgradSt);
			//print("s",l,s,x,y,"testhere",xtry,ytry);
			//Overlay.drawRect(xtry, ytry, 2,2);
			//Overlay.show;
			   v=getPixel(xtry, ytry);
			 if (v==0) setPixel(xtry, ytry, 8); // Avoid zapping outline pixels in line 322
			  //print("v",v);
			   if (v==2) {  // start at an end of line point
			   		xBegin=xtry;
			   		yBegin=ytry;
			   		fnd=1;
			   		//print("     found end 2               ",l,"  ",xBegin,yBegin);
			   		//break;
			   } // if v==2
			   if (fnd==1) break;
			   	//setPixel(xtry, ytry, 8);
			} // for y loop
		   Overlay.show;
		if (fnd==1) break;
			//setPixel(xtry, ytry, 8);
			} // x loop
			if (fnd==1) break; // have start point
		} // s loop   search for start point	
		//print("have found end",fnd,"pos",xBegin,yBegin);

// Now Follow Sequential pixels along membrane and measure
run("Select None");
changeValues(8, 8, 0); // remove markers from search for start
run("Duplicate...", " ");
rename("BinaryPerim" + l);
UseforPerimImage=getTitle();
//changeValues(8, 255, 0);

setMinAndMax(0, 100); // display range
changeValues(1, 255, 255); // make binary    
//jjjjj

// how many binary pixels in line- need to create arrays of correct size
// option could make a larger array and then trim afterwards
wd=getWidth();// size
ht=getHeight();
getRawStatistics(nPixels, mean);
nLine=1+round(mean*wd*ht/255); // +1 to accommodate start pos
	//print("nLine",nLine);
	PerimXarray=newArray(nLine); // empty array, 0s
	PerimYarray=newArray(nLine);
		lastNotes=lastNotes + 1;
	atnline=lastNotes;
	  //bbbbbbbbbbbb
	setResult("Notes", atnline, d2s(round(nLine),0));// N of pixels in PM
	//setResult("Notes", atnline, round(nLine));// N of pixels in PM
	row=0; // was row=0, but added mean and std to results
		PerimXarray[row]=xBegin;// start position  
	   	PerimYarray[row]=yBegin;
	   		if(l==0) setResult("xposO", row+2, xBegin);// start position
			if(l==0) setResult("yposO", row+2, yBegin); // start position
				setPixel(xBegin,yBegin, 64);// zap current - prevent stepping backwards
	row=1;// was row=0, but added mean and std to results
	
// Have Start Point and Number of Locations -FIND COORDINATES   
// create arrays to hold values
//     or read GP values
AtX=xBegin;// start at - location Need- single pixel
AtY=yBegin;
   //print("start",0,AtX,AtY);
listpos=1;
// search around point, changes to coordinates
Xarray=newArray(1, -1, 0, 0, 1, -1, -1,  1);// 8  test pos - orthog first then diagonal
Yarray=newArray(0 , 0, 1,-1, 1, -1,  1, -1);
setBatchMode(1);
selectImage(UseforPerimImage);
for (r=1;r<=nLine;r++) { // search for connected binary pixels    -find line
fndPixel=0; // flag  check that a pixel on the line is actually found
for (g=0;g<8;g++) { // search 8 neighboring pixel
			//print(r,"   g",g);
	selectImage(UseforPerimImage);// not needed    check
	   run("Remove Overlay");
	// check orthogonal - prefer orthogonal move
	xTst=AtX + Xarray[g];
	yTst=AtY + Yarray[g];
	v=getPixel(xTst,yTst);// read pixel value
			//print("        g",v);
			//print(r,"   g",g,v,xTst,yTst);
	if(v>100) { // part of line
		fndPixel=1; // found pixel on line=xTst;
			//print(r,"   g",g,v,xTst,yTst,fndPixel);
		AtX=xTst;
		AtY=yTst;
		setPixel(AtX,AtY, 10+g);// zap current - prevent stepping backwards
		Overlay.drawRect(AtX,AtY, 0.5,0.5);
		 //print("new",r,row,AtX,AtY);
		 	PerimXarray[row]=AtX;// save location in array
	   		PerimYarray[row]=AtY;
	   		if(l==0) setResult("xposO", row+2, AtX);// add pos to Results
			if(l==0) setResult("yposO", row+2, AtY);// add pos to Results
			row=row+1; 	
		break; // leave for g loop  look around
	} // if
	//if (fndPixel==0) print("    line pixel not found around",r, AtX,AtY);
	//if (fndPixel==1) break; // found next binary pixel, escape g loop
	//if (fndPixel==0) print("    line pixel not found around",r, xTst,yTst);
 } // for g loop    follow line, looking at neighboring pixels
 // if (fndPixel==0) print("    line pixel not found around",r, AtX,AtY);

  if (fndPixel==0)	break; // end of perimeter   escape r loop
} // for r loop
Overlay.show;

  if (fndPixel==0) { // run out of binary pixels,  why?
	// reduce size of the two coordinate arrays - delete 0,0 location
	PerimXarray=Array.slice(PerimXarray,0,row);
	PerimYarray=Array.slice(PerimYarray,0,row);
  }
// print("pixels in line",l,nLine,"actual",row);

setBatchMode(0);

/* print out line coordinates for one line  l==0 central
print("points on line ");
for (rd=0;rd<PerimXarray.length;rd=rd+1) {
	if (l==0) print(l,rd,PerimXarray[rd],PerimYarray[rd]);
}
*/

// Have sequential points, now read channels and add to Results
setBatchMode(1);
//selectImage(ImageMeasureName);
selectImage(title);

setSlice(1); 
print(" ");
Array442=newArray(PerimXarray.length);
for (rd=0;rd<PerimXarray.length;rd=rd+1) {
	v=getPixel(PerimXarray[rd],PerimYarray[rd]);
		//xmark=PerimXarray[rd];
		//ymark=PerimYarray[rd];
		//Overlay.drawRect(xmark,ymark, 1, 1);
	v=v-BKgrnd442;// subtract background
	//print(rd,PerimXarray[rd],PerimYarray[rd],v);
		Array442[rd] = v;
	setResult(columnNameSlice1, rd+2, d2s(v,0));
	 //print(rd,PerimXarray[rd],PerimYarray[rd],v);
	//if(l==-3) setResult("xpos", rd, xmark);// check on positions
	//if(l==-3) setResult("ypos", rd, ymark);
} // for rd
Overlay.show;

setSlice(2);
Array483=newArray(PerimXarray.length);
for (rd=0;rd<PerimXarray.length;rd=rd+1) {
	v=getPixel(PerimXarray[rd],PerimYarray[rd]);
	v=v-BKgrnd483;// subtract background
		Array483[rd] = v;
	setResult(columnNameSlice2, rd+2, d2s(v,0));
} // for rd 

// GP-values
print(" ");
GPlabel="GP_" + l;
GParray=newArray(PerimXarray.length);
for (rd=0;rd<PerimXarray.length;rd=rd+1) {
	GP=(Array442[rd] - Array483[rd]) / (Array483[rd] + Array442[rd]);
	//print("rd",rd,"GP",GP,Array483[rd],Array442[rd]);
	setResult(GPlabel, rd+2, d2s(GP,3));
	GParray[rd]=GP;
} // for rd 
Array.getStatistics(GParray, GPmin, GPmax, GPmean, GPstdDev);
//print(GPlabel,GParray.length, GPmin, GPmax, GPmean, GPstdDev);
setResult(GPlabel, 0,  d2s(GPmean,3));
setResult(GPlabel, 1,  d2s(GPstdDev,3));

//setBatchMode(0);
setBatchMode("exit and display");
//close();
} // L lOOP vary distance from membrane....................................................

setForegroundColor(255,255,255);

// add note to Results table 
additionalNote=getString("Note in Result Table", "good");
setResult("Notes", 12,"  "+additionalNote); // user's comments
updateResults();
selectWindow("Results");

// Option to save Results & ROIs  - in original folder
SaveYesNo=getBoolean("Save Results & ROIs");
// save in   UseDirectory
if(SaveYesNo==1) {
	SaveText=getString("ROI's new name", ROItitle);
	//print("SaveText",SaveText);	
	SaveBasic=replace(title, ".tif", "");// 
	SaveBasic=replace(SaveBasic,"Images","Results");
		//print(title,"  ",SaveBasic);
	
} // save
