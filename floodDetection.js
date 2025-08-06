//***********************New and improved biweekly flood maps*********************************//

// Defining the AOI
var districts = ee.FeatureCollection("FAO/GAUL/2015/level2");

//assam
// var aoi = districts.filter(ee.Filter.and(
//   ee.Filter.eq("ADM1_NAME", "Assam"),
//   ee.Filter.or(
//     ee.Filter.eq("ADM2_NAME", "Kamrup"),
//     ee.Filter.eq("ADM2_NAME", "Barpeta"),
//     ee.Filter.eq("ADM2_NAME", "Nalbari")
//   )
// ));


//Kerala
var region = districts.filter(ee.Filter.and(
  ee.Filter.eq("ADM1_NAME", "Kerala"),
    ee.Filter.or(
    ee.Filter.eq("ADM2_NAME", "Ernakulam"),
    ee.Filter.eq("ADM2_NAME", "Kottayam"),
    ee.Filter.eq("ADM2_NAME", "Alappuzha")
  )
));

/////////////////////////////////////
//created a new defination for AOI(required due to large size of raster)
//defining the region bounds in meters
var regionBoundsSize = 6000;
var aoi;

// Parameters for shifting
var shiftX_meters = -11000; // negative = left (west), positive = right (east)
var shiftY_meters = -10000;  // positive = up (north), negative = down (south)

if (regionBoundsSize) {
  var centroid = region.geometry().centroid();

  // Convert meters to degrees (approximate)
  var metersToDegrees = function(meters) {
    return ee.Number(meters).divide(111320);  // Approx. at equator
  };

  var shiftX_deg = metersToDegrees(shiftX_meters);
  var shiftY_deg = metersToDegrees(shiftY_meters);

  var coords = centroid.coordinates();
  var lon = ee.Number(coords.get(0));
  var lat = ee.Number(coords.get(1));

  var shiftedLon = lon.add(shiftX_deg);
  var shiftedLat = lat.add(shiftY_deg);

  var shiftedCentroid = ee.Geometry.Point([shiftedLon, shiftedLat]);

  // Buffer and bound to make square AOI
  aoi = shiftedCentroid.buffer(regionBoundsSize / 2).bounds();
} else {
  aoi = region.geometry().bounds();
  print('Using default AOI bounds for export.');
}

Map.centerObject(region, 10);
////////////////////////////////////////////////////////////


//*****************display the aoi boudary********************
  var empty = ee.Image().byte();
  var aoiOutline = empty.paint({
    featureCollection: ee.FeatureCollection(aoi),
    color: 1,
    width: 3 // Defines the line width of the outline.
  });
  
  Map.addLayer(aoiOutline, {palette: 'yellow'}, 'Export Bounding Box');
//************************************

//CONFIGS
var years = ee.List.sequence(2018, 2024);
var monsoonStart = 6, monsoonEnd = 9;
var threshold = -16;



/**
 * Returns Image Collection of water mask over the entire year(monsoon season)
 * uses a threshold to mask out water and non water
 * if there are more than 1 images in 2-weeks then it takes the mean of them
 * if there is no image for the biweek(rare) then a complete 0 image is taken(TODO - we can take the prev available image or mean of the prev and next image)
 * Param - year
 **/
function biWeeklyMasks(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, monsoonStart, 1);
  var end   = ee.Date.fromYMD(year, monsoonEnd, 30);
  
  var nBiWks = end.difference(start, 'week').divide(2).floor();

  // Create a list of start dates for each bi-weekly period
  var starts = ee.List.sequence(0, nBiWks.subtract(1))
    .map(function(i) {
      return start.advance(ee.Number(i).multiply(2), 'week');
    });

  // Collect the images and return them
  var biWeeklyMasksCollection  =  ee.ImageCollection(starts.map(function(w0) {
    w0 = ee.Date(w0);
    // The end of the period is 2 weeks after the start
    var w1 = w0.advance(2, 'week');

    // Get Sentinel‑1 VV data for the 2-week period
    var col = ee.ImageCollection('COPERNICUS/S1_GRD')
      .filterBounds(aoi)
      .filterDate(w0, w1)
      .filter(ee.Filter.eq('instrumentMode','IW'))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
      .select('VV');

    var isCollectionEmpty = col.size().eq(0);
    var biWeekly = ee.Image(ee.Algorithms.If(isCollectionEmpty, ee.Image().select(), col.mean()));
    var mask = biWeekly.lt(threshold).unmask(0);
    var biWkNum = w0.difference(start, 'week').divide(2).floor();

    return mask
      .rename("water")
      .set({
        'year': year,
        'biweek': biWkNum,
        'system:time_start': w0.millis(),
        'has_data': ee.Algorithms.If(isCollectionEmpty, 0, 1)
      })
      .clip(aoi);
  }));
  
  return biWeeklyMasksCollection;
}
var biWeeklyMasksByYear = years.map(biWeeklyMasks); // calling the above function and storing masks for all the years


//*********************
//*********************
//*********************
//*********************

/**
***************DATA IS COLLECTED AT THIS POINT. STARTING THE CLASSIFICATION***********************
**/


/**
 * STEP 1...
 * Classifying into perennial, non water and unclassified pixel.
 * This classification is done using the entire data collected and can be updated using the CONFIG variable
 * There will be a further classification post selecting the year and week.
 * That classification will further classify the unclassified pixels into perennial,flood,seasonal for the given year and week
 * 
 **/
function classifyPerennialAndNonWater(years, monsoonStart, monsoonEnd, aoi, threshold, perennialThreshold) {

  //Calculate the total number of bi-weekly periods across all years.
  var totalBiWeeks = years.map(function(year) {
    var start = ee.Date.fromYMD(year, monsoonStart, 1);
    var end = ee.Date.fromYMD(year, monsoonEnd, 30);
    return end.difference(start, 'week').divide(2).floor();
  }).reduce(ee.Reducer.sum());

  //Sum the bi-weekly water masks to get water presence count.
  var waterPresenceSum = ee.Image(0);
  for (var i = 0; i < years.size().getInfo(); i++) {
    var yearCollection = ee.ImageCollection(biWeeklyMasksByYear.get(i));
    if (yearCollection.size().getInfo() > 0) {
      waterPresenceSum = waterPresenceSum.add(yearCollection.reduce(ee.Reducer.sum()).unmask());
    } else {
      waterPresenceSum = waterPresenceSum.add(ee.Image(0).toByte());
    }
  }

  // proportion of time water was present
  var waterPresenceProportion = waterPresenceSum.divide(ee.Image.constant(totalBiWeeks));
  
  // Classifying the pixels
  var perennialWater = waterPresenceProportion.gte(perennialThreshold).rename('perennial_water'); //perennial threshold is 0.85
  var nonWater = waterPresenceProportion.eq(0).rename('non_water'); // if the water presence was 0 across all images then its a permanent non water pixel

  var classification = ee.Image(0)
    .where(perennialWater, 1)
    .where(nonWater, 2);

  return classification.clip(aoi);
}

//storing the classification 
//This is the initial classification. Further classification requires selecting 
//the year and date.
var perennialAndNonWaterClassification = classifyPerennialAndNonWater(years, monsoonStart, monsoonEnd, aoi, threshold, 0.85);

//!This can be removed later(for testing purposes)!
Map.addLayer(perennialAndNonWaterClassification, { palette: ['000000', '0000FF','00BB00'], min:0, max:2 }, 'Perennial and Non-Water');





/**
*STEP 2...
* creating a function for classifying the unclassified pixels
* This function will be triggered when the inputs are provided and generates a fully classified flood map for 
* the provided period
* 
**/
function displaySelectedMask() {
  var weekFreq = 0.6; // This is now bi-weekly frequency
  var yearFreq = 0.8;

  var selectedYear = yearSelect.getValue();
  var selectedBiWeek = ee.Number(biWeekSelect.getValue()).toInt();


  var selectedYearCollection = ee.ImageCollection(biWeeklyMasksByYear.get(years.indexOf(selectedYear)));

  // Filter the ImageCollection for the selected bi-week
  var selectedMask = selectedYearCollection
    .filterMetadata('biweek', 'equals', selectedBiWeek)
    .first();
    
  if (selectedMask) {
    // Since we confirmed selectedMask exists, we cast it to an ee.Image
    selectedMask = ee.Image(selectedMask);
    var selectedBiWeekStartDate = ee.Date.fromYMD(selectedYear, monsoonStart, 1).advance(selectedBiWeek.multiply(2), 'week');
    print('Processing data for bi-week starting:', selectedBiWeekStartDate.format('YYYY-MM-dd'));


    var monsoonBiWeeks = years.map(function(year) {
      var start = ee.Date.fromYMD(year, monsoonStart, 1);
      var end = ee.Date.fromYMD(year, monsoonEnd, 30);
      return end.difference(start, 'week').divide(2).floor();
    });
    var totalMonsoonBiWeeks = ee.Number(monsoonBiWeeks.reduce(ee.Reducer.sum()));

    // Find water frequency for the "selected bi-week" across all years
    var waterFrequencyThisBiWeek = ee.ImageCollection(years.map(function(year) {
      var yearCollection = ee.ImageCollection(biWeeklyMasksByYear.get(years.indexOf(year)));
      var biWeekImage = yearCollection.filterMetadata('biweek', 'equals', selectedBiWeek).first();
      return ee.Algorithms.If(biWeekImage, ee.Image(biWeekImage), ee.Image(0).selfMask());
    })).reduce(ee.Reducer.sum()).divide(years.size());

    // Find water frequency for the "selected year" across all bi-weeks
    var waterFrequencyThisYear = selectedYearCollection.reduce(ee.Reducer.sum()).divide(totalMonsoonBiWeeks);

    //non water pixel for the provided week
    var weeklyNonWaterMask = selectedMask.eq(0).rename('weekly_non_water'); //the water is present on the pixel for the given bi-week
    var unclassifiedPixels = perennialAndNonWaterClassification.eq(0);
    var newNonWaterPixels = unclassifiedPixels.and(weeklyNonWaterMask);

    var alreadyClassified = perennialAndNonWaterClassification;
    alreadyClassified = alreadyClassified.where(newNonWaterPixels, 2);

    var seasonalWater = waterFrequencyThisBiWeek.gte(weekFreq).and(waterFrequencyThisYear.lt(yearFreq)).rename('seasonal_water').and(alreadyClassified.eq(0));
    var floodWater = waterFrequencyThisBiWeek.lt(weekFreq).and(waterFrequencyThisYear.lt(yearFreq)).rename('flood_water').and(alreadyClassified.eq(0).and(seasonalWater.not()));
    var newPerrenialWater = waterFrequencyThisYear.gte(yearFreq).rename('new_seasonal_water').and(alreadyClassified.eq(0).and(seasonalWater.not().and(floodWater.not())));

    var finalClassification = alreadyClassified.where(seasonalWater, 3).where(floodWater, 4).where(newPerrenialWater, 5);
    
    //****** additional filters *******//
    var floodWaterVector = createFloodWaterVectors(finalClassification, aoi, 200000,weekFreq,yearFreq);
    
    
    //rendering the image
    Map.layers().set(1, ui.Map.Layer(finalClassification.clip(aoi), {
      palette: ['000000', '0000FF', '00BB00', 'yellow', 'red', '72A1ED'],
      min: 0,
      max: 5
    }, 'Water Classification'));

    // Export.image.toDrive({
    //   image: finalClassification,
    //   description: 'finalClassification_biweekly',
    //   fileNamePrefix: 'finalClassification_biweekly',
    //   region: aoi.geometry(),
    //   scale: 30,
    //   crs: 'EPSG:4326',
    //   maxPixels: 1e13
    // });

  } else {
    // This block runs if no image was found for the selection.
    print('No data found for the selected year and bi-week.');
    Map.layers().set(1, ui.Map.Layer(ee.Image().rename('Water Classification'), {}, 'Water Classification', false));
  }
}


//***************DEFINING THE HELPER FUNCTIONS **********************************

/**
 * Converts the class of pixel into vector polygons
 * then filters out the polygons based on  min area
 * Params - classificationImage (the final classified image)
 *        - floodClassValue (the class to vectorize)
 *
 **/
function createFloodWaterVectors(classificationImage, aoi, minAreaSqm,weekFreq,yearFreq) {
  
  Map.clear(); // clear the prev image
  Map.centerObject(aoi, 5);
  
  var exportScale = 30;
  print('Native classification scale (m):', classificationImage.projection().nominalScale());
  print('Using export scale (m):', exportScale);
  

  // Step 1: Isolate flood pixels (class 3 or 4).
  var floodMask = classificationImage.eq(3).or(classificationImage.eq(4)).selfMask();
  
  
  // Step 2: Smudge/dilate the mask to generalize the area and connect nearby pixels.
  var smudgedMask = floodMask.focal_max({
    radius: 30,
    units: 'meters'
  }).selfMask();

  // Step 3: Vectorize the mask to create polygons from the raster data.
  var floodVectors = smudgedMask.reduceToVectors({
    geometry: aoi,
    scale: 30,
    geometryType: 'polygon',
    labelProperty: 'class',
    maxPixels: 1e13
  });
  
  // --- DEBUGGING STEP ---
  // Check if any vectors were created before area calculation.
  print('Number of polygons BEFORE area filtering:', floodVectors.size());

  // Step 4: Compute the area of each polygon and store it as a property.
  floodVectors = floodVectors.map(function(feature) {
    var area = feature.geometry().area({maxError: 1});
    return feature.set({'area_m2': area});
  });

  // Step 5: Filter out polygons that are smaller than the specified minimum area.
  var filteredVectors = floodVectors.filter(ee.Filter.gte('area_m2', minAreaSqm));
  
  // --- DEBUGGING STEP ---
  // Check if any vectors remain after filtering.
  print('Number of polygons AFTER area filtering:', filteredVectors.size());
  
  // Render the filtered vector polygons to the map for visualization.
  Map.addLayer(filteredVectors, {color: 'red'}, 'Filtered Flood Polygons');

  //creating a new vector after filtering the polygons(It would be difficult to compare polygon)
  var emptyImage = classificationImage.multiply(0).byte();

  // "Paint" the polygons onto this correctly-projected empty image.
  // Pixels covered by a polygon get a value of 1. This avoids reprojection issues.
  var finalFloodRaster = emptyImage.paint({
    featureCollection: filteredVectors,
    color: 1 // This is the value assigned to the pixels within the polygons.
  }).clip(aoi);

  // Render the final rasterized flood map. The .selfMask() will hide 0 values.
  Map.addLayer(finalFloodRaster.selfMask(), {palette: ['#0000FF']}, 'Final Flood Raster (Masked)');
  
  // Add the raster WITHOUT the mask to confirm it's being created correctly.
  // You should now see blue areas where your polygons are.
  Map.addLayer(finalFloodRaster, {min: 0, max: 1, palette: ['black', 'blue']}, 'Final Flood Raster (No Mask)');


  // Step 7: Export the final FeatureCollection to your Google Drive.
  var yearStr = String(yearFreq).replace(/\./g, '_');
  var weekStr = String(weekFreq).replace(/\./g, '_');
  
  var rasterFileName = 'flood_raster_' + yearStr + '_week' + weekStr;
  var rasterDescription = 'Flood_Raster_' + yearStr + '_week' + weekStr;

  
  
  
  Export.image.toDrive({
    image: finalFloodRaster,
    description: rasterDescription,
    folder: 'keshav_sparsh/FloodMaps',
    fileNamePrefix: rasterFileName,
    region: aoi,
    scale: exportScale,
    crs: classificationImage.projection()
  });


  // Step 8: Return the final FeatureCollection.
  return filteredVectors;
}

//***************************
//***************************
//***************************
//***************************
//***************************

// -------------------- Start of the UI part ----------------------------
var panel = ui.Panel({ style: { width: '300px' } });
var yearSelect = ui.Select({
  items: years.getInfo().map(function(year) {
    return { label: String(year), value: year };
  }),
  placeholder: 'Select Year'
});

// The number of bi-weekly periods is about half the number of weeks.
// June 1 to Sept 30 is ~17 weeks, so there are ~8 bi-weekly periods.
var biWeekSelect = ui.Select({
  items: ee.List.sequence(0, 8).getInfo().map(function(week) {
    return { label: 'Bi-Week ' + String(week), value: week };
  }),
  placeholder: 'Select Bi-Week'
});
var displayButton = ui.Button({ label: 'Show Bi-Weekly Water' });
var weekLayer = ui.Map.Layer(ee.Image().rename('Bi-Weekly Water'), { palette: ['000000', '0000FF'] }, 'Bi-Weekly Water');
Map.layers().push(weekLayer);

displayButton.onClick(displaySelectedMask);

panel.add(ui.Label('Select Year and Bi-Week to Display'));
panel.add(yearSelect);
panel.add(biWeekSelect);
panel.add(displayButton);
ui.root.add(panel);

// --- Legend ---
var legendPanel = ui.Panel({
  style: { position: 'bottom-left', padding: '8px 15px' }
});
var legendTitle = ui.Label({
  value: 'Water Classification Legend',
  style: { fontWeight: 'bold', fontSize: '16px', margin: '0 0 4px 0' }
});
legendPanel.add(legendTitle);
var palette = [ '0000FF', '00BB00', 'FFFF00', 'FF0000'];
var names = ['Perennial Water', 'Perennial Non Water/Non-Water', 'Seasonal Water', 'Flood Water'];
for (var i = 0; i < palette.length; i++) {
  var colorBox = ui.Label({
    style: { backgroundColor: '#' + palette[i], padding: '8px', margin: '0 0 4px 0' }
  });
  var description = ui.Label({ value: names[i], style: { margin: '0 0 4px 6px' } });
  legendPanel.add(ui.Panel({
    widgets: [colorBox, description],
    layout: ui.Panel.Layout.Flow('horizontal')
  }));
}
panel.add(legendPanel);
