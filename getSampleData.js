// ****************************************************************************************
// Description: Extract EVI2 time series from Landsat 5, 7, 8, and 9 products for 
// point locations. This is to supplement the blsp package in R
// Author: J. Gao, Ian McGregor
// Last updated: June 2022
// ****************************************************************************************

// STEP 1: DEFINE VARIABLES AND YOUR DESIRED POINTS
// *************** Define Variables ************************
// Please adjust this as you need. Note L9 has no data prior to 2021
var startDate = "1984-01-01";
var endDate = "2022-06-01";

var sensor = 
  ["LANDSAT/LT05/C02/T1_L2",    //L5 Collection 2 SR
    "LANDSAT/LE07/C02/T1_L2",   //L7 Collection 2 SR
    "LANDSAT/LC08/C02/T1_L2",   //L8 Collection 2 SR
    "LANDSAT/LC09/C02/T1_L2"    //L9 Collection 2 SR
  ];
  
// columns for the output csv
var bandNames = ['satellite', 'date', 'lon', 'lat', 'id', 'evi2', 'QA_PIXEL']

// Labels for exporting the csv to google drive
var taskDescription = "exportEVI2" //description of task in GEE
var folder = "Bayesian_LSP" //folder to export to, single string
var fileName = "sampleData" //name of file, e.g. this will be sampleData.csv

// *************** Define Variables ************************
// site locations can be either be specified by coordinates...

// A sample list of coordinates
var points = ee.List([
  [-71.700975, 43.945733], // HQ
  [-71.725558, 43.952022], // 1B.
  [-71.728739, 43.950917], // 4B.
  [-71.731661, 43.958947], // 4T.
  [-71.731831, 43.949128], // 5B.
  [-71.736869, 43.957400], // 5T.
  [-71.742178, 43.955547], // 6T.
  [-71.765647, 43.928928], // 7B.
  [-71.769733, 43.918358]  // 7T.
  ]);

// Create feature collection from points
var samp_pts = ee.FeatureCollection(points.map(function(p){
  var point = ee.Feature(
    ee.Geometry.Point(p), {
      id: points.indexOf(p)
    }
  );
  return point;
}));

// ...Or, upload a csv file with columns "id", "longitude", "latitude" columns, 
// and add it to this code as an asset
//var samp_pts = ee.FeatureCollection(table);

print("Sample Points", samp_pts);

// *************** Define Individual Functions ************************
/**
* Function to mask clouds based on the pixel_qa band of Landsat data.
* @param {ee.Image} image Input Landsat SR image
* @return {ee.Image} Cloudmasked Landsat SR image
*/
function mask_landsat_sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 4);
  var cloudsBitMask = (1 << 3);
  // Get the pixel QA band.
  var qa = image.select('QA_PIXEL');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

/**
* Function to calculate EVI2.
* @param {nirBand} string NIR band
* @param {redBand} string Red band
* @param {image} a Landsat image
* @return {ee.Image} Landsat SR image with EVI2 band
*/
function calcEVI2(nirBand, redBand){
  var calc = function(image){
    var nir = image.select([nirBand]).multiply(0.0000275).add(-0.2);
    var red = image.select([redBand]).multiply(0.0000275).add(-0.2);
    var evi2 = image.expression(
                '2.5 * ((nir - red) / (1 + nir + 2.4 * red))', {
                  'nir': nir,
                  'red': red
                });
    return image.addBands(evi2.rename("evi2"));
  };
  return(calc);
}

/**
* Wrapper function for all steps of the code: for each sensor, load data,
* filter by date and points, cloudmask, get evi2, then identify wanted values
* @param {sensorName} a Landsat image with evi2 band
* @return {ee.FeatureCollection} a collection of features
*/
function wrapper(sensorName){
  // Make sure the name is a string
  var sens = ee.String(sensorName);
  
  // Identify which Landsat it is (matters for NIR and Red band selection below)
  var landsatType = sens.match('[8-9]');

  var evi2Bands = ee.Algorithms.If(landsatType.size().eq(0), 
                  ee.List(["SR_B3", "SR_B4", "QA_PIXEL"]),
                  ee.List(["SR_B4", "SR_B5", "QA_PIXEL"]));
  var nir = ee.List(evi2Bands).get(0);
  var red = ee.List(evi2Bands).get(1);

  // Create the image collection, then filter by date and points, then 
  // apply cloudmask and evi2 calculation
  var IC = ee.ImageCollection(sens);

  var fullTS = IC.filterDate(startDate, endDate)
              .filterBounds(samp_pts)
              .select(evi2Bands)
              .map(mask_landsat_sr)
              .map(calcEVI2(nir, red));
  
  // Get the wanted values for each point for each image
  function getValsEachImage(img){
    
  // Now we get the wanted values per point. Specifically, we make
  // a feature for each point, whose properties are the values we
  // want to export. This returns a feature collection
    var getVals = samp_pts.map(function(pt){
      var ID = pt.get('id');
      var thispoint = pt.geometry();
      var thiscoord = thispoint.coordinates();
  
      var date = img.date().format('yyyy-MM-dd');
      var value = img.select(['evi2', 'QA_PIXEL'])
                  .reduceRegion(ee.Reducer.first(), thispoint);
    
      var output = ee.Feature(thispoint,
                  {date: date, satellite: sens, id: ID,
                    lat: thiscoord.get(0), lon: thiscoord.get(1)
                  })
                  .set(value);
    
      return(output);
    });
    return(getVals);
  }

  var out = fullTS.map(getValsEachImage);
  return(ee.FeatureCollection(out.flatten()));
}

// *************** Run the main code ************************
var loop = sensor.map(wrapper);
print("Main output", loop);

// Here, loop is a list because 'sensor' is a list, so we need to merge
var combined = loop[0].merge(loop[1]).merge(loop[2]).merge(loop[3]);
print("First element example", combined.first());

// Export to csv
Export.table.toDrive(
  {collection: combined,
    description: taskDescription,
    folder: folder,
    fileNamePrefix: fileName,
    selectors: bandNames
  }
);
