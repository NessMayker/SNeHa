# SNeHa

This readme will detail the notebook order and processes used in this project. 

0.
We first locate the available H-alpha files using both MUSE and ESO-MPG and download and clean the OSC catalog.

	Notebook 0.MUSEIFU_Ha_Images:
	This notebook opens the MUSE files and locates the SNe within them, plotting them in a folder.
	This notebook also writes the MUSE galaxies and their filenames in ../Data/MUSEdata.csv
	The MUSE files have this address: /data/fourier/sun.1608/PHANGS/MUSE/DR2.1/MUSEDAP/copt/

	Notebook 0.ESO-MPG_Ha_Images:
	This notebook opens the ESO-MPG files and locates the SNe within them, plotting them in a folder.
	This notebook also writes the ESO-MPG galaxies and their filenames in ../Data/NBdata.csv
	The ESO-MPG files have this address: /home/mayker.1/Desktop/SNeHaLargeData/HaSUB_wcomb_corr/

	Notebook 0.OpenSupernovaCatalogCleanerCode:
	This notebook takes the original downloaded Open Supernova Catalog csv file and "cleans" it by 		removing SNe without discovery dates, and averaging the reported RAs and Decs.
	Then converts the desired entries into a region file.

1.
We then create a database for the project.

	Notebook 1.ProjectTable:
	This notebook pulls the MUSE and ESO-MPG data tables and joins them, compiling a database for this 		project.


2.
We then find the recent supernovae that occur in both the MUSE and ESO footprints. We plot these and take rudimentary intensity measurements. !we still need to measure errors and record beamsize for ESO files!

	Notebook 2.SNe Mapper:
	This notebook reads in the cleaned OSC file and searches through ESO & MUSE map images, looking for 		OSC objects that have gone off in their footprint. 
	It then converts the ra & dec coordinates to map pixels and measures the pixel intensity there.
	This notebook prints the OSC Object data on the '2.SNeHaMasterCat.txt'. 
	This notebook then turns the 2.SNeHaMasterCat.txt file into a dataframe, where the non SN objects 		are removed and then the cleaned table is printed to '2.SNeHaMasterCatClean.csv'.
	This cleaned table is used to plot all of the galaxy images and their respective SNe 




