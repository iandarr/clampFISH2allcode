

An example pipeline is described to get from an experiment's
raw data & user-input specification tables --> extractedData --> plots/tables & example images

In this example, the short experiment identifier is "E159" (with the longer identifier "E159_ampclick_and_scan") and these plots, tables, and example image are specifically from 'part 1' of this experiment (relating to amplification and +/-click). The plots to be generated will demonstrate amplification and click (in Fig 1, and also in the supplement).

--------------------------------------------- Example folder structure ---------------------------------------------
The folder structure for the files relevant to this example is as follows, divided at the highest level between "rawData" directory and the "paper" directory. Since the rawData directory has the image files, it is a very large size. As such, the subfolders of the rawData directory are often on external hard disks, while the paperDirectory is often on the C: drive. These absolute paths will need to be specified at the start of the scripts. For example:
rawDataDirectory='/Volumes/IAND_04/'
paperDirectory='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper/'

/../../rawData
	E159_ampclick_and_scan
		part1_8well_ampclick
			... / subfolders /
				XYdata.csv 
				dapi001.tif
				gfp001.tif
				gfpa001.tif
				gfpb001.tif
				gfpc001.tif
				gfpd001.tif
				gfpe001.tif
				tmr001.tif				
				tmra001.tif				
				tmrb001.tif				
				tmrc001.tif				
				tmrd001.tif				
				alexa001.tif				
				alexaa001.tif				
				alexab001.tif				
				alexac001.tif				
				alexad001.tif				
				cy001.tif				
				cya001.tif				
				cyb001.tif				
				cyc001.tif				
				cyd001.tif				
				trans001.tif				
				dapi002.tif
				... (more .tif files)
				data001.mat
				data002.mat
				... (more data00X.mat files)
				mergedImgs
					Merged001.tif
					Merged001_cp_masks.png
					Merged002.tif
					Merged002_cp_masks.png
					... (more files)
/../../paper
	experiments
		E159part1_amp_click_Tcond.xlsx
		E159part1_amp_click_ImagingDetails.xlsx
	scripts
		rawdataprocessingScripts
			E159part1_rawdataprocessing.m
		extractionScripts
			E159part1_extract.m
		plotScripts
			E159part1_plot.m
		showImagesScripts
			E159_exampleImages_amp.m
		functions
	extractedData
		E159_ampclick_scan
			part1_ampclick
				... (some intermediate Tspots and Tcell tables, some in subfolders) ...
				Tspots4AllNormalized.csv
				Tcell2All.csv
				Tcomb.csv
				TcombConsise.csv
				
	plots
		amplification_and_click
			(a plot input file would be here, but none exists for this example)
			E159_amplification_click_UBC.eps
			E159_amplification_4genes.eps
			E159_click_4genes.eps
	exampleImages
		amplification
			amplification_exampleImagesInputs.xlsx
			amplification_exampleImagesFinal.xlsx
			R1_10X.tif
			R1_20X.tif
			R1_60X.tif
			R2_10X.tif
			... (more images)
	figures
		

		
--------------------------------------------- Example workflow ---------------------------------------------
The overall workflow looks like this:
rawData (initial ND2 files and corresponding Tiffs) ---rawdataprocessingScript-->  rawData (possibly renamed Tiffs, +/- Merged YFP+DAPI tiffs, +/- png masks, and data00X.mat rajlabimagetools files)
rawData ---extractionScript (initial code section)--->  extractedData (initial Tspots.csv and Tcells.csv
extractedData  ---extractionScript (later code section)--->  extractedData


- Use E159part1_rawdataprocessing.m to:
	rename tif files to rajlabimagetools using '... _ImagingDetails.xlsx' file
	merge DAPI00X.tif and YFP00X.tif files into an RGB image, output to the mergedImgs subfolder with files named Merged00X.tif
	use cellpose to segment Merged00X.tif files, with the outputs named Merged00X_cp_masks.png
	generate rajlabimagetools data00X.mat based on segmentations files. (Note: In experiments with manual segmentation, the previous two 'merged00X.tif' and cellpose steps are not done). The 60X segmentations are used to make segmentations for the corresponding 10X and 20X data.
	process data00X.mat files (aTrous filter to find regional maxima, set thresholds, exclude slices). Each of these steps results in modified data00X.mat files.
- Use E159part1_extract.m to:
	extract spot data from the data00X.mat files and put them into tables in the extractedData directory.
- Use E159part1_plot.m to:
	plot the results, which takes its data from tables from the extractedData folder
	