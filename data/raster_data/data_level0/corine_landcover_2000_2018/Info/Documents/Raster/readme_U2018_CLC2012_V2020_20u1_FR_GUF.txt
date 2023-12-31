﻿README - CLC2012 France - French Guyana
----------------------------------------

Layer name: U2018_CLC2012_V2020_20u1_FR_GUF

Version 2020 20u1 (dated 02/2020): First Update of Version 20 

Database covers 1 country. See https://land.copernicus.eu/user-corner/technical-library/clc-country-coverage-1990-2018-v20u1 for full information about the coverage of this version.

Vector version is created in ESRI GeoDatabase format. 
Some redundant lines between neighbouring polygons with the same code are still present in database, but only as result of persisting ‘adaptive tilling’ procedure limitation of ArcGIS.
Polygons under MMU 25 Ha can be present along national borders and along 'adaptive tilling' tiles boundaries. Polygons (slivers) under 0.1 Ha can be present only along 'adaptive tilling' tiles boundaries.

Raster version is created in 100m resolution in GeoTiff format, using CELL_CENTER method for the rasterizing. 

File naming coventions:

Filename is composed of combination of information about update campaign, data theme and reference year and version specification (including release year and release number) 
UPDATE CAMPAIGN | THEME | REFERENCE YEAR | RELEASE YEAR | RELEASE NUMBER 

UPDATE CAMPAIGN:
•	update campaign refers to the year of the CLC campaign, when the file content was last modified or updated 
THEME:
•	theme refers to specific CLC layer type: clc - refers to status layer, cha - refers to change layer
REFERENCE YEAR: 
•	reference year refers to the year for which CLC information included in the file were mapped
RELEASE YEAR:
•	release year refers to the year of CLC data release 
RELEASE NUMBER:
•	release number refers to a sequential numbering of CLC data releases. Initially it is named as a beta version (e.g. 20b2) until the new complete final version is ready covering the full territory (e.g. 20). Any subsequent minor update is recorded as incremental number in suffix (e.g. 20u1)


Examples:
U<update campaign>_<theme><reference year>_V<release year>_<release number>

status file name: U2006_CLC2000_V2020_20b2' means
U2006 - last modification / update of the file is from the 2006 mapping campaign 
CLC2000 - file contains CLC status data for reference year 2000
V2020_20b2 - file was released in 2020 as second beta version (incomplete) of version 20


change file name: file name 'U2018_CHA1218_V2020_20u1' means
U2018 - last modification / update of the file is from the 2018 mapping campaign 
CHA1218 - file contains CLC change data between reference years 2012 and 2018 
V2020_20u1 - file was released in 2020 as a first update of final version 20

Change layers in raster
----------------------
In order to keep 8bit resolution for raster version, change layers are divided into two files. 

Vector file name:  
U2018_CHA1218_V2020_20u1

Raster files name: 
U2018_CHA1218_12_V2020_20u1 - file contains consumption part of change (‘from’ code)
U2018_CHA1218_18_V2020_20u1 - file contains formation part of change (‘to’ code)

Overseas countries and territories - Overseas departments (DOMs)
----------------------------------------------------------------
Overseas countries and territories are delivered as separate package to keep consistency between vector and raster delivery. Structure of the delivery and file naming conventions are the same except DOM specification (country code and DOM code) suffix.
UPDATE CAMPAIGN | THEME | REFERENCE YEAR | RELEASE YEAR | RELEASE NUMBER | COUNTRY CODE | DOM CODE

Example:
U<update campaign>_<theme><reference year>_V<release year>_<release number>_<country code>_<DOM code>

file name: U2006_CLC2000_V2018_20b2_FR_GLP' means
U2006 - last modification / update of the file is from the 2006 mapping campaign 
CLC2000 - file contains CLC status data for reference year 2000
V2018_20b2 - file was released in 2018 as second beta version (incomplete) of version 20
FR - file is for French territory 
GLP – DOM name is Guadeloupe
 

In order to keep 8bit resolution for change layers, they are divided into two files. 
U2018_CHA1218_12_V2020_20u1 - file contains consumption part of change (‘from’ code)
U2018_CHA1218_18_V2020_20u1 - file contains formation part of change (‘to’ code)

Change log (from previous release):

Version 20u1
------------ 
Release date: 24-02-2020
Main purpose of the release: Maintenance / Correction of final CLC2018 data.
Changes from previous release (20): 
•	File naming conventions simplified and better described. New file naming convention has been introduced based on user feedback on version 20. Filename is composed of combination of information about update campaign, data theme and reference year and version specification (including release year and release number). See https://land.copernicus.eu/user-corner/technical-library/clc-file-naming-conventions-guide-v20u1 
•	The French DOMs are provided in separate databases (files both for vector and raster version of data).
•	All raster layers are back in 8 bit GeoTIFF. Modification is introduced based on the user feedback on version 20. In order to keep 8 bit resolution for raster change layers, they are divided into two files - representing consumption (from) and formation (to) part of change.
See https://land.copernicus.eu/user-corner/technical-library/clc-country-coverage-1990-2018-v20u1 for full information about the coverage of this version.

Version 20
----------
Release date:1-05-2019
No changes for CLC2006 FR DOMs layers
See https://land.copernicus.eu/user-corner/technical-library/clc-country-coverage-1990-2018-v20 for full information about the coverage of this version.

Version 20b2
------------
Release date: 19-12-2018
Coding for rasters changed - rasters contains directly CLC codes/changes values. So status rasters produced as 16-bit depth, change rasters produced as 32-bit depth 
See https://land.copernicus.eu/user-corner/technical-library/clc-country-coverage-v20b2/view for full for full information about this version coverage

Version 20b1 
------------
Release date 16-11-2018
Only raster file re-rasterized by CELL CENTRE method was delivered.

Prepared by
GISAT 02/2020
For more information contact tomas.soukup@gisat.cz, miroslav.kopecky@gisat.cz