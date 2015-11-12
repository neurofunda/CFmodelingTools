%--------------------------------------------%
To attach the cc-pRF menus to the mrVISTA GUI:

1. Add the following line to 'mrLoadRet/Analysis/analysisMenu.m', just below the call to the retinotopic model submenu:
ccmAnalysisMenu(analysismenu,vw);

2. Add the following line to 'mrLoadRet/file/filesMenu.m', just below the call to the retinotopic model submenu:
ccmLoadMenu(fileMenu, view);

3. Add the following line to 'mrLoadRet/Plots/plotMenu.m', just below the call to the retinotopic model submenu:
ccmPlotMenu(h, view);


%--------------------------------------------%
To run the connective field analysis, first select the appropriate datatype and load at least two rois. 
You might then want to restrict these rois to layer-1.
It is then possible to start the analysis via the analysis menu: Analysis >> Cortico-Cortical Model >> Run 

The analysis will first ask you to select the source rois (usually V1-Layer1) and the target roi. 
You can only select one target roi; if you want to run the analysis for multiple rois, you should first combine them.

The analysis will then start with computing a huge matrix that contains all pairwise distances between voxels in the
source roi. It will use this matrix to generate the connective field models and then the model predictions. 

Sometimes, I get an out-of-memory error. This is usually due to the fact that the arrays containing the time-series
are too large. Converting them into single precision may help. If you're working under Windows, turning off the 3GB 
switch could also help; this allows matlab to allocate more continous memory slots.

The core of the analysis uses matlab's lscov function, which takes most of the analysis' time. If you are using a 
machine with multiple cores, lscov will use all cores as long as the matrix contains more than 40K elements. The lscov
function also doesn't seem to care much about how many target voxels we wish to analyse; the length of the time-series
seems to be much more important.

Once the analysis finished, it writes out a ccModel*.mat file containing the results of the analysis. You can load this
file via File >> Cortico-Cortical Model >> Select and Load Cortico-Cortical Model. This will load the connective field 
size in the amplitude view and the variance explained in the coherence view. Optionally, you can then derive the pRF
eccentricity and polar-angle values via the loaded connective field model via File >> Cortico-Cortical Model >> 
Load Derived Maps. You will then be asked to select a pRF model from which the eccentricity and polar-angle values in
the source region will be used to derive the same parameter values for the voxels in the target region. Make sure, 
therefore, to have the appropriate rois loaded.

With the connective field model loaded, it is possible to plot several things. The plot options can be found in the plots menu
under Cortico-Cortical Model.

