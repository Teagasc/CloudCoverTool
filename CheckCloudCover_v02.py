# -*- coding: utf-8 -*-
"""
Spyder Editor

The aim of this script is to assess the amount of useful (cloud free) satellite
imagery on a parcel/field base

It can either take centorids or full polygons as an input for the parcels and
will read Sentinel-2 A and B images. 

"""

# =============================================================================
# Import required packaged
# =============================================================================

import os, shutil, glob, argparse, sys, tempfile, tarfile, zipfile
import numpy as np
import pandas as pd
import rasterio as rio
import geopandas as gpd
from rasterio.mask import mask
from shapely.geometry import mapping


# =============================================================================
# Define functions
# =============================================================================

class params():
    
    '''
    A class containing the necessary parameters for the tool.
    '''
    
    platform = None
    FeatureDirectory = None
    identifier = None
    ImageDirectory = None
    
    def __init__(self, platform = None, FeatureDirectory = None, 
                 identifier = None, ImageDirectory = None, parser = None):  
        if parser:
            
            self.readParser(parser)
        else:
            self.platform = platform
            self.FeatureDirectory = FeatureDirectory
            self.identifier = identifier
            self.ImageDirectory = ImageDirectory
        
    def readParser(self, parser):
        self.platform = parser.platform[0]
        self.FeatureDirectory = parser.Features[0]
        self.identifier = parser.index[0]
        self.ImageDirectory = parser.ImageDirectory[0]

def any_in(a, b): 
    '''
    A handy little function that will check two lists for coninciding items
    Returns True/False
    '''
    return any(i in b for i in a)
    
def extract_by_feature(inRaster, extractGDF, featureID, *args):
    '''
    The function extracts raster values based on an input geopandas 
    geodataframe. 
    The funciton uses the rasterio.mask function to crop the input raster with
    each element of the extraction geodataframe. Additional arguments can be 
    used to exclude extracts that contain the specified values 
    
    It has been tested with points and polygons.
    '''
    
    # Create an empty dictionary populate with the extracted arrays
    outArrays = dict()     
    
    # Using rasterio to read the raster
    with rio.open(inRaster) as rioRaster:
        # Reproject the geodata grame to the coordinate reference system of the 
        # raster
        projGDF = extractGDF.to_crs(rioRaster.crs.to_dict())
        
        # Iterate over all elements of the geodataframe
        for i in range(len(projGDF)):
            # Convert geodataframe element to geoJSON to be readible by 
            # rasterio mask function
            geojson = [mapping(projGDF.geometry[i])]
            try:
                # Use the rasterio mask function to create the raster extract
                # and the raster transformation. The nodata value is used as a
                # filler for final array. 
                out_image, out_transform = mask(rioRaster, geojson, crop=True, 
                                                nodata = 25)
                
                # exclude the extract if all cells have a no data value
                if np.unique(out_image) != 25:
                    if args:
                        # If additional arguments are given, exclude tiles that 
                        # include the given values
                        if any_in(args, out_image):
                            continue   
                        else:
                            pass
                        
                    outArrays[extractGDF[featureID][i]] = out_image
            except:
                continue
            
    if len(outArrays) == 0:
        return None
    else:        
        return outArrays
    
def extract_tar(inTarFile, outFolder, file = None):
    '''
    Uses the tarfile package to extract .tar and .tgz files. 
    Can be used to extract single files  
    '''
    if not os.path.exists(outFolder):
        print('Output directory did not exist, created {} in {}'.format(os.path.basename(inTarFile), os.path.dirname(inTarFile)))
        os.makedirs(outFolder)

    if not os.path.isfile(inTarFile):
        raise(ValueError(".tar file does not exists"))
        
    tar = tarfile.open(inTarFile)
    if file:
        tar.extract(file, outFolder)
    else:
        tar.extractall(path=outFolder)
    tar.close()
    
def extract_zip(inZipFile, outFolder, file = None):
    '''
    Uses the zipfile package to extract .zip archives files. 
    Can be used to extract single files  
    '''
    
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
        
    if not os.path.isfile(inZipFile):
        raise(".zip file does not exists")
    
    zip_ref = zipfile.ZipFile(inZipFile,"r")    
    if not file:
        files = [inZip for inZip in zip_ref.namelist() if '.gfs' not in inZip]        
        zip_ref.extractall(outFolder, files)
    else:
        zip_ref.extract(member = file, path = outFolder)
    zip_ref.close()
    
def createMaskfromBits(qaFile, b1, b2):
    '''
    Reads bit-encoded rasters and decodes them into masks representing 
    specified bit(s). This can be used to extract cloud or other masks from 
    the quality layers provided with Landsat and Venus products.
    
    The values in the main script are set to the standard encoding of the 
    products. Please consult with product manuals to identify the specific 
    values encoded.
    '''
    qaBand = qaFile.read()
    width_int = int((b1 - b2 + 1) * "1", 2) 
    return ((qaBand >> b2) & width_int).astype('uint8')

def writeArrayToRaster(inArray, outRaster, crs, transform, driver = "GTiff"):
    '''
    A function to write a numpy array to a GeoTIFF raster. Requires a 
    coordinate reference system and a transfomation provided.
    '''
    with rio.open(outRaster, 'w', driver = driver,
                            height = inArray.shape[1], width = inArray.shape[2],
                            count = 1, dtype=str(inArray.dtype),
                            crs = crs,
                            transform = transform) as rasterFile:
        rasterFile.write(inArray)
    
    return rasterFile

def readSentinel2QualityMasks(s2Folder, s2mask = 'MSK_CLOUDS_B00'):
    '''
    Sentinel-2 Quality Masks are provided as .gml files. This function reads 
    user defined masks and reads them as geopandas geodataframes
    '''
    maskFile = glob.glob(os.path.join(s2Folder, "GRANULE", os.listdir(os.path.join(s2Folder, "GRANULE"))[0], "QI_DATA", "*{}*.*".format(s2mask)))[0]
    
    # geopandas.read_file cannot extract the crs from .gml files. Instead the
    # crs is taken from the pre-view image jpeg200 
    crsFile = glob.glob(os.path.join(s2Folder, "GRANULE", os.listdir(os.path.join(s2Folder, "GRANULE"))[0], "QI_DATA", '*.jp2'))[0]
    
    # Use rasterio to extract crs from pre-view image
    with rio.open(crsFile) as crs_layer:
        crs = crs_layer.crs.to_epsg()  
    # read mask to geopandas dataframe and assign crs from pre-view image
    maskData = gpd.read_file(maskFile, driver = 'gml')
    maskData.crs = crs
    
    # return the mask as a unary union geodataframe. The reference system is
    # also returned for later use
    return(maskData.unary_union, maskData.crs)

def read_FC_to_GDF(fullFCPath):
    '''
    Currently not in use. 
    Will be used to read features from ESRI file geodatabases
    '''
    
    # Split full path to feature class into geodatabase/feature dataset path, 
    # and feature class name
    inputGDB = os.path.dirname(fullFCPath)
    inputFC = os.path.basename(fullFCPath)
    
    outGDF = gpd.read_file(inputGDB, layer = inputFC)
    return(outGDF)
        
def readFeaturesToGDF(inFeatures):
    '''
    Reads shapefiles (either all shapefiles in a given folder, or a single 
    file) and converts them to a geopandas geodataframe.
    '''
    
    if os.path.isdir(inFeatures):
        # if a directory is provided, the tool will read all shapefiles in it
        inData = [os.path.join(FeatureDirectory, inShape) for inShape in os.listdir(FeatureDirectory) if inShape.endswith('.shp')]
        featuresList = list()
        # all shapefiles in the directory are appended to a new list 
        # as geodataframes
        for inFile in inData:
            featuresList.append(gpd.read_file(inFile).explode())
        # Using the pandas concat function the list is concatenated into a 
        # single geodataframe
        featuresGDF = pd.concat(featuresList, ignore_index = True)
        # The index is reset to avoid duplications
        featuresGDF = featuresGDF.reset_index()
    else:
        # if a single file is provided it is directly read to a geodataframe
        featuresGDF = gpd.read_file(inFeatures)
        
    if featuresGDF.crs == None:
        raise(Exception("Input feature class does not have a spatial reference system defined."))
    return(featuresGDF)
        
def checkS2CloudCover(dataFolder, inFeatures, identifier):
    '''
    Checks if a set of features is obscured by cloud cover on a Sentinel-2 
    image. Cloud cover is obtained either from the Sentinel-2 QI files or by 
    using the fmask algorithm. 
    '''
    
    try:
        featureOverlay = dict()
        cloudMask, crs = readSentinel2QualityMasks(dataFolder, s2mask = 'MSK_CLOUDS_B00')
    except ValueError:                               
        sensorFootprint, crs = readSentinel2QualityMasks(dataFolder, s2mask = 'MSK_DETFOO_B01')
        
        inFeatures = inFeatures.to_crs(epsg = crs)  
        for idx, point in inFeatures.iterrows():
            if point.geometry.within(sensorFootprint):
                featureOverlay[point[identifier]] = os.path.basename(dataFolder)
        return(featureOverlay) 
    
    sensorFootprint, crs = readSentinel2QualityMasks(dataFolder, s2mask = 'MSK_DETFOO_B01')
    
    inFeatures = inFeatures.to_crs(epsg = crs) 
    for idx, point in inFeatures.iterrows():
        if point.geometry.within(sensorFootprint) and not point.geometry.within(cloudMask):
           featureOverlay[point[identifier]] = os.path.basename(dataFolder)

    if len(featureOverlay) > 0:
        return(featureOverlay)    
    else: return(None)
    
def checkVenusCloudCover(dataFolder, inFeatures):
    '''
    Checks if a set of features is obscured by cloud cover on a Venus image. 
    Cloud cover is obtained from the Venus QI layer 
    ''' 
    cloudLayer = os.path.join(dataFolder, [folder for folder in os.listdir(dataFolder) if folder.endswith('CLD.DBL.TIF')][0])
    with rio.open(cloudLayer) as tempLayer:    
        cloudMaskArray = createMaskfromBits(tempLayer, 0, 0)
        crs = tempLayer.crs.to_dict() 
        
        transform = tempLayer.transform
    cloudMask = os.path.join(tempDirectory, "tempCCMask.tif")
    writeArrayToRaster(cloudMaskArray, cloudMask, crs, transform)        
    featureOverlay = extract_by_feature(cloudMask, inFeatures, identifier, -999, 1) 
    return(featureOverlay)
    
def LandsatCloudCover(dataFolder, inFeatures):
    cloudLayer = glob.glob(dataFolder + '\\*qa*.*')[0]
    band = 5 if 'pixel_qa' in cloudLayer else 4
    with rio.open(cloudLayer) as tempLayer:    
        cloudMaskArray = createMaskfromBits(tempLayer, band, band)
        crs = tempLayer.crs.to_dict() 
        transform = tempLayer.transform
    cloudMask = os.path.join(tempDirectory, "tempCCMask.tif")
    writeArrayToRaster(cloudMaskArray, cloudMask, crs, transform)  
    featureOverlay = extract_by_feature(cloudMask, inFeatures, identifier, 1) 
    return(featureOverlay)
    
def main():
    finalFeatureList = dict(site = list(), image = list())
    inFeatures = readFeaturesToGDF(FeatureDirectory)
    
    imageIdentifier = {'Sentinel2': 'S2',
                       'Venus': 'VENUS',
                       'Landsat8': 'LC8'}
    
    imageArchive = [os.path.join(ImageDirectory, folder) for folder in os.listdir(ImageDirectory) if imageIdentifier[platform] in folder]

    for image in imageArchive:
        print("Checking {}...".format(os.path.basename(image)))
        if os.path.isdir(image):
            if len(os.listdir(image)) > 0:
                DECOMPRESSED = False
                pass
            
        elif image.endswith('.tar') or image.endswith('.tgz'):
            extract_tar(image, tempDirectory)
            DECOMPRESSED = True
            
        elif image.endswith('.zip'):
            extract_zip(image, tempDirectory)
            DECOMPRESSED = True
            
        if platform == 'Venus':
            if DECOMPRESSED: 
                dataFolder = os.path.join(tempDirectory, [folder for folder in os.listdir(tempDirectory) if folder.endswith('.DIR')][0])
            else: dataFolder = image
            featureOverlay = checkVenusCloudCover(dataFolder, inFeatures)
            
        elif platform == "Landsat8":
            if DECOMPRESSED: 
                dataFolder = tempDirectory
            featureOverlay = checkS2CloudCover(dataFolder, inFeatures, identifier = identifier)
        
        elif platform == "Sentinel2":
            if DECOMPRESSED: 
                dataFolder = os.path.join(tempDirectory, [folder for folder in os.listdir(tempDirectory) if folder.endswith('.SAFE')][0])
                
            else: dataFolder = image
            featureOverlay = checkS2CloudCover(dataFolder, inFeatures, identifier = identifier)
                   
        if featureOverlay:
            for key, value in featureOverlay.items():
                finalFeatureList['site'].append(key)  
                finalFeatureList['image'].append(os.path.basename(image))  
         
        shutil.rmtree(tempDirectory)
                
    finalFeatureFrame = pd.DataFrame(data = finalFeatureList)  
    finalFeatureFrame.to_csv(os.path.join(ImageDirectory, "CloudFreeFeatures_{}.csv".format(platform)))#
    
# =============================================================================
#       Setting parameters
# =============================================================================

if __name__ == "__main__":
    
    sys.path.append('c:\programdata\miniconda\envs\scripts')

    parser = argparse.ArgumentParser(description='Check cloud-cover of multiple satellite images at specific locations defined by features')
    parser.add_argument('Features', type = str, nargs = 1,
                    help='Shapefile or directory containing a set of shapefiles')
    parser.add_argument('ImageDirectory', type = str, nargs = 1,
                    help='Directory containing satellite images that should be tested for cloud cover')
    parser.add_argument('platform', type = str, nargs = 1, choices = ['Sentinel2', 'Venus', 'Landsat'],
                        help = 'Select satellite platform (important for selecting the right cloud cover)')
    parser.add_argument('-i', '--index', type = str, nargs = 1, metavar = '', default = 'index',
                        help = 'In case of multiple features is is advised to provide a unique feature identifier field to identify the feature that was recognised  as cloud free. Will default to a generic "index"')
    parser.add_argument('-f', '--fmask', type = bool, nargs = 1, default = False, metavar = '',
                        help = 'Calcluate Sentinel-2 cloud mask using the fmask algorithm. Default: False ')

# =============================================================================
#     Parsing input arguments
# =============================================================================
 
    args = parser.parse_args()

    config = params(parser = args)
        
    platform = config.platform
    FeatureDirectory = config.FeatureDirectory
    identifier = config.identifier
    ImageDirectory = config.ImageDirectory
    tempDirectory = tempfile.mkdtemp(dir = ImageDirectory, prefix = 'CloudCover')
    print(FeatureDirectory)

# =============================================================================
# Execute main function
# =============================================================================

    main()
                   
