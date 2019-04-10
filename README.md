# CloudCoverTool

The tool will scan an archive of satellite imagery (Sentinel-2 Level-1C; Landsat 8 Level-1C, SR; Venµs SR) 
and return all images that according to the product specific cloud mask are cloud free at a user defined
area of interest.

## Introduction
The tool takes a single feature, a feature class, or a directory containing multiple feature classes 
(any format readable by geoandas) as a user defined area of interest. It will then craws through a user defined 
satellite image archive and extract the cloud mask supplied with the products. 
Each cloud mask will be overlaid with the area(s) of interest and all features not intersecting with cloud cover
will be extracted. 
The tool output is a csv file containing a feature identifier (either user defined or a simple index) and a list
of images that were identified as showing no cloud cover at that specific location.

## Script
The tool is written in Python and can be used from a command line.
Inputs are:<br/>
<b> Required </b><br/>
Features defining the area of interest<br/>
Location of the imagery archive<br/>
Satellite Platform (Sentinel-2, Landsat-8, Venus)<br/>
<b> Optional </b><br/>
-i Feature identifier<br/>
