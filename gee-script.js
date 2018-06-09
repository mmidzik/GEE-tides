//Harmful Algal Bloom Analysis from MODIS  //Maya Midzik: mmidzik@gmail.com 
//FES 754 
 
/* This  code differentiates ride tides/harmful algal blooms (HABs)  from other near coastal waters and phytoplankton blooms 
using  Moderate Resolution Imaging Spectroradiometer (MODIS)  satellite imagery and the methods of Hu et al. (2013) 
 
MODIS first must be processed from Level 1a MODIS data 
(raw radiance counts) to Level 1b files (geolocated  
at-aperture radiances) using NASA's SeaWiFS  Data  Analysis    System (SeaDAS). These files were then atmospherically corrected  using the SeaDAS l2gen code, generating remote sensing reflectance,  chl-a, particulate backscattering  using  the  QAA 
and normalized water-leaving radiance. 
 
Inputs to the HAB MODIS model are: 
1.Enhanced RGB (ERGB) image 
2. Chl-a 
3.Fluorescence Line Height (FLH) 
4.Particulate  backscattering coefficient  ratio
calculated from: bbp,551/bbp,Morel*/   
//--------Raster inputs from SeaDAS------- 
 
//1. Enhanced RGB  (ERGB) image bands: 
//Normalized water-leaving radiance (nLw) at: 547, 488 and 443 nm 
//var MODIS_ERGB = ee.Image('GME/images/ 14493593916934751216-11613121305523030954').select('b1', 'b2', 'b3');  var MODIS_ERGB = ee.Image('GME/images/ 14493593916934751216-17488941626782682984').select('b1', 'b2', 'b3');   
//2. Chl-a bands: 
var MODIS_CHL = ee.Image('GME/images/ 14493593916934751216-04594044960761723980').select('b1', 'b2') 
 
//3. FLH bands (to be calculated) 
var MODIS_FLH = ee.Image('GME/images/ 14493593916934751216-11853667273131550346').select('b1','b2','b3') 
 
 
 //4. 
//--a) particulate backscatterig coefficient at 551nm using QAA 
var MODIS_bbp = ee.Image('GME/images/ 14493593916934751216-08551536342581976716').select('b1') 
//--b) bands for bbp Morel calculation 
var MODIS_CHl = ee.Image('GME/images/ 14493593916934751216-16143158689603361093').select('b1') 
 
//-------Processing input rasters—— 
 
//-- Center and zoom to region of interest based on individual scene   //Scene based centering coded to allow application to any initial inputs  //Define scene extent of images 
var SceneExtent = MODIS_CHL.select('b1'); 
var feature = ee.Feature(SceneExtent) 
var GetBounds = feature.bounds() 
//Get centroid of scene extent 
var center = GetBounds.centroid(null, 'EPSG:4326') 
var latlong = center.geometry() 
//Convert centroid to numerical coorinates 
var x= (latlong.coordinates().get(0).getInfo()) 
var y= (latlong.coordinates().get(1).getInfo()) 
var latlongnum= ee.List(latlong.coordinates()) 
print(latlongnum) 
Map.setCenter(x,y,7) 
 
//---Add land/water mask to entire image to allow better visualization of   //model input layers 
var landwatermask = ee.Image('MODIS/MOD44W/ MOD44W_005_2000_02_24').select('water_mask');  print(landwatermask) 
var vizParams = { 
min: 0, 
max: 1, 
palette: ['808080','000000'], 
}; 
Map.addLayer(landwatermask, vizParams, 'Land water Mask') 
 
//1. --- Visualize Enhanced RGB Image 
var imageERGB = MODIS_ERGB.visualize({bands: ['b1', 'b2', 'b3'], max: 20});  Map.addLayer(imageERGB, null, 'ERGB') 
 
 //2. --- Calculate Chlorophyll from MODIS OC3 equation  //NASA OC3: 
//Rrs1 = blue wavelength Rrs (e.g., 443, 490, or 510-nm)  //Rrs2 = green wavelength Rrs (e.g., 547, 555, or 565-nm)  //X = log10(Rrs1 / Rrs2) 
//chlor_a = 10^(a0 + a1*X + a2*X^2 + a3*X^3 + a4*X^4) 
var ReflRatio = MODIS_CHL.select('b1').divide(MODIS_CHL.select('b2'))  var logRefl = ReflRatio.log10() 
var Chl_OC3 = MODIS_CHL.expression( 
'10**(0.2500 + -2.4752*R + 1.4061*(R**2) + -2.8233*(R**3) + 0.5405*(R**4))', {  'R':logRefl 
}) 
var ChlViz = {min: 0, max: 10, palette: ['0000FF', '00FF00', 'FF0000']}; 
//Visualize CHL layer 
Map.addLayer(Chl_OC3, ChlViz, 'CHL_A') 
 
//3. --- Calculate FLH from respective loaded wavelengths 
//First caluclate band ratio of respective bands 
var BandRatio=ee.Number((746.3-676.7)/(746.3-665.1)) 
var FLH_baseline= MODIS_FLH.select('b1').add((MODIS_FLH.select('b2').subtract(MODIS_FLH.select('b1'))).multipl y(BandRatio)) 
//Calculate FLH for the entire scene 
var FLH_scene = MODIS_FLH.select('b3').subtract(FLH_baseline); 
//Visualize the FLH layer 
var FLHViz = {min: 0, max: 1, palette: ['0000FF', '00FF00', 'FF0000']};  Map.addLayer(FLH_scene, FLHViz, 'FLH') 
 
 
//4. --- Calculate bbp Morel from loaded SeaDAS bands 
//First, mask out any Chl levels under 1.5 for calculation to exclude nonproductive waters 
var ChlMask = MODIS_CHl.gt(1.5) 
var Masked_chl = MODIS_CHl.mask(ChlMask) 
 
//A. -- Add bbp from original SeaDAS file (bbp using Quasi-Analytical alogrithm) 
var Masked_bbp = MODIS_bbp.mask(ChlMask) 
 
//B. -- Calculate bbpMorel using the Morel (1988) algorithm 
var bbpMorel = MODIS_CHL.expression( 
'(0.3*powChl)*(0.002+(0.02*(0.5-0.25*logChl)))', {  'powChl': Masked_chl.pow(0.62), 
'logChl': Masked_chl.log10() 
})   

 //C. -- Calculate bbp ratio of bbpMorel/bbpQAA 
var BbpRatio = Masked_bbp.divide(bbpMorel) 
 
//Mask out values of zero (pixels without full ratio calculation)  var RatioMask = BbpRatio.gt(0) 
var BbpRatioMasked = BbpRatio.mask(RatioMask) 
//visualize bbpratio layer 
var ratioViz = {min: .1, max: 2, palette: ['FF0000', '0000FF','008000' ]};  Map.addLayer(BbpRatioMasked, ratioViz, 'ratio') 
 
//----Use given parameters to define classification regions of ocean color 
 
//1. --- Define dark and bright ERGB areas for classification 
//define the 25th percentile of reflectance for Red (448) and Green (557) bands   
var ERGB_Red = imageERGB.select('vis-red').reduceRegion({ 
reducer: ee.Reducer.percentile([25], ['25th'], null, null, null),  geometry: feature.geometry(), 
scale: 30, 
maxPixels: 1e9 
}); 
var ERGB_Green = imageERGB.select('vis-green').reduceRegion({ 
reducer: ee.Reducer.percentile([25], ['25th'], null, null, null),  geometry: feature.geometry(), 
scale: 30, 
maxPixels: 1e9 
}); 
 
var Rednumber = ERGB_Red.get('vis-red'); 
var Greennumber= ERGB_Green.get('vis-green'); 
 
//Use the percentiles to extract the high and low reflectance of ERGB bands  var ERGBhighRefl = function (number_ith,Initial_scene) { 
var ERGB_ith_number = ee.Number(number_ith) 
var ERGB_ith_mask= Initial_scene.gt(ERGB_ith_number) 
var ERGB_ith_masked = Initial_scene.mask(ERGB_ith_mask)  return ERGB_ith_masked 

var ERGBlowRefl = function (number_ith,Initial_scene) { 
var ERGB_ith_number = ee.Number(number_ith) 
var ERGB_ith_mask= Initial_scene.lt(ERGB_ith_number) 
var ERGB_ith_masked = Initial_scene.mask(ERGB_ith_mask)  return ERGB_ith_masked 
} 

  
//Apply to red (448) and green (557) bands 
var RedPercentileupperRefl= ERGBhighRefl(Rednumber,(imageERGB.select('vis-red')))  var RedPercentilelowerRefl= ERGBlowRefl(Rednumber,(imageERGB.select('vis-red')))  var GreenPercentileupperRefl= ERGBhighRefl(Greennumber,(imageERGB.select('vis- green'))) 
var GreenPercentilelowerRefl= ERGBlowRefl(Greennumber,(imageERGB.select('vis- green'))) 
 
//Define bright and dark ERGB areas for final evaluation 
var ERGB_bright=RedPercentileupperRefl.and(GreenPercentileupperRefl) 
var ERGB_dark=RedPercentilelowerRefl.and(GreenPercentilelowerRefl) 
 
//2. ---No additional processing of Chl scenes 
 
//3. --- Define areas of high and low FLH values 
//In order to find areas of high FLH, first define the 90th percentile of calculated FLH  
var FLH_75th =FLH_scene.reduceRegion({ 
reducer: ee.Reducer.percentile([90], ['75th'], null, null, null),  geometry: feature.geometry(), 
scale: 30, 
maxPixels: 1e9 
}); 
print(FLH_75th) 
var number75th = FLH_75th.get('b3'); 
 
 
///Final.-----Defining categories of water and creating Chl. error    
/*In order to combine water classificaitons into a single 
values of each individual image must be set to a binary grid.  (ie. masked masked values must be converted into 0's in oder  for bands to be combined) */ 
 
var imagemask = function (Initial_image) { 
var imagemask=Initial_image.eq(1) 
var maskedimage=Initial_image.mask(imagemask)  return maskedimage 
} 
 
var blank= ee.Image(0) 
var Wherefunction= function (Initial_image){ 
var output=blank.where(Initial_image.eq(1),1) 
return output  } 
     

 ///---Compute chlorophyll error matrix (CDOM-rich water miss-classified  //as Chl-a from NASA OC3 algorithm) 
var scaledFLH = FLH_scene.multiply(10) 
var ChlError = MODIS_CHl.subtract(scaledFLH)  Map.addLayer(ChlError, null , 'Chlorophyll Error') 
 
//----Give categorical definitions of each water classification 
//based on previously defined parameters, and apply above 
//functions to eliminate masked values 
var SedimentRich = Wherefunction(imagemask(ERGB_bright.and(FLH_scene.gt(ee.Number(number75th)))))  var ShallowClear = Wherefunction(imagemask(ERGB_bright.and(FLH_scene.lt(ee.Number(number75th)))))  var CDOMRich = Wherefunction(imagemask(ERGB_dark.and(FLH_scene.lt(ee.Number(number75th)))))  var PhytoRich = Wherefunction(imagemask(ERGB_dark.and(FLH_scene.gt(ee.Number(number75th)))))  var HABPoor = Wherefunction(imagemask(PhytoRich.and(BbpRatio.gt(.8)))) 
var HABRich = Wherefunction(imagemask(PhytoRich.and(BbpRatio.lt(.8)))); 
 
 
//Merge all categories into a single 'classification' Image 
//This image will have 5 bands, each with binary values 
var mosaic=ee.Image.cat([SedimentRich,ShallowClear,CDOMRich,HABPoor,HABRich])  print(mosaic, 'mosaic') 
 
//Covert binary layer into single band layer with values 1-5 
//corresponding to class definitions 
var CatImage = mosaic.expression( 
'Class1 + 2*Class2 + 3*Class3 + 4*Class4 + 5*Class5', {  'Class1': mosaic.select('constant'), 
'Class2': mosaic.select('constant_1'), 
'Class3': mosaic.select('constant_2'), 
'Class4': mosaic.select('constant_3'), 
'Class5': mosaic.select('constant_4'),  }); 
var igppalette= {min: 0, max: 5, palette: ['FFFFFF','00FFFF','0000FF','00FF00', 'FF0000']}   
//Mask out all zero values of final image to reveal Land/Wayer Mask 
var CatImageMask = CatImage.gt(0) 
var FinalClassified = CatImage.mask(CatImageMask) 
//Visualize final classification image  Map.addLayer(FinalClassified, igppalette, 'Final Classified Water')