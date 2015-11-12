#CFmodelingTools

In the paper, [Connective field modeling](http://www.ncbi.nlm.nih.gov/pubmed/23110879), Koen Haak introduces a novel method to analyze cortico-cortical functional connectivity between retinotopic cortices. Here I show how to implement the code with mrVista to compute connective field (CF) models in matlab 2011b.

The code provided in this repository is intended to retrieve already computed pRF model parameters, ROIs times series and compute CF models. The idea is to use the function `pRFandCFdata` to collect, compute and store, in an easy to handle structure, all relevent parameters. The function depends on the following inputs:

A loaded mrVista session (view) with loaded source and target ROIs.
Path to a folder containing the pRF model.
Datatype numbers (This is vestigial for now. It can be useful to compute many dataTypes at once).
Variance explained threshold for plotting pRF maps (optional).
Number of frames and TR.
The output of this function is the structure `pRFandCFdata_dataType_sourceROI_targetROI`, which contains:

* Source and target ROI time series.
* Associated pRF variance explained (VE).
* Associated pRF size (Sigma).
* Associated Visual field coordinates (Eccentricity and Polar angle).
* Associated Cortical coordinates.
* Cortical distances between all pair of nodes (voxels).

And the sub-structure "cf", which contains the CF parameters:

* Variance explained.
* Corrected variance explained.
* Sigma.
* Rss, Raw Rss, beta, npoints, ntrends.
* CF center coordinates.
* Cartesian and polar visuotopic coordinates.
* Target voxel coordinates.
* Connection index (each target ROI index and its associated CF index in the source ROI).


This function also computes the power spectrum of the detrended time series, plots them and save the images as pdf files. The power spectrum is computed using Chronux.

See example inside `pRFandCFdata`.

nicolas.gravel@gmail.com 11-11-2015 University of Groningen, The Netherlands.

Always cite your sources!

References:

[Haak KV, Winawer J, Harvey BM, Renken R, Dumoulin SO, Wandell Ba, et al. Connective field modeling. NeuroImage. 2013;66:376–384.](http://www.ncbi.nlm.nih.gov/pubmed/23110879)

[Gravel NA, Harvey B, Nordhjem B, Haak KV, Dumoulin SO, Renken R, et al. Cortical connective field estimates from resting state fMRI activity. Frontiers in Neuroscience. 2014;8(October):1–10.](http://www.ncbi.nlm.nih.gov/pubmed/25400541)(