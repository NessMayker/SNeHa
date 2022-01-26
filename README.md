# SNeHa

This readme will detail the notebook order and processes used in this project. 

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
	This notebook takes the original downloaded Open Supernova Catalog csv file and "cleans" it by removing SNe 		without discovery dates, and averaging the reported RAs and Decs.
	Then converts the desired entries into a region file.

We then create a database for the project.

	Notebook 1.ProjectTable:
	This notebook pulls the MUSE and ESO-MPG data tables and joins them, compiling a database for this project.



