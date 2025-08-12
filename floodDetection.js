// ==================== CONFIGURATION PARAMETERS ====================
var CONFIG = {


  // Time parameters
  years: [2018, 2019, 2020, 2021, 2022, 2023, 2024],
  monsoonStart: 6,
  monsoonEnd: 9,
  
  // Classification thresholds
  threshold: -16, //sentinel 1 VV threshold
  perennialThreshold: 0.85,
  weekFreq: 0.6,
  yearFreq: 0.8,
  
  // Processing parameters
  minAreaSqm: 200000,
  regionBoundsSize: 0, // Set to 0 to use full AOI bounds
  
  // Export parameters
  exportToAsset: false,
  exportToDrive: true,
  assetPath: "users/your_username/flood_maps/",
  driveFolder: "keshav_sparsh/FloodMaps",
  
  // Batch processing mode
  batchMode: false, // Set to true for automated batch processing
  batchYears: [2023, 2024], // Years to process in batch mode
  batchBiWeeks: [3, 4, 5] // Bi-weeks to process in batch mode
};






/////////////////************STARTING THE CODE *****************/////////////////////////////
// ==================== REGION SETUP ====================
function setupRegion() {
  
  var districts = ee.FeatureCollection("FAO/GAUL/2015/level2");
  var region = districts.filter(ee.Filter.and(
    ee.Filter.eq("ADM1_NAME", "Kerala"),
      ee.Filter.or(
      ee.Filter.eq("ADM2_NAME", "Ernakulam"),
      ee.Filter.eq("ADM2_NAME", "Kottayam"),
      ee.Filter.eq("ADM2_NAME", "Alappuzha"),
      ee.Filter.eq("ADM2_NAME", "Pattanamtitta")
    )
  ));


  var aoi;
  if (CONFIG.regionBoundsSize) {//clip the boundary based on the clipping box
    aoi = region.geometry().centroid().buffer(CONFIG.regionBoundsSize / 2).bounds();
    print('Using centered export box of size (m):', CONFIG.regionBoundsSize);
  } else {
    aoi = region; //return the entire region
    print('Using default AOI bounds for export.');
  }
  
  return aoi;
}


// Initialize AOI
var aoi = setupRegion();
Map.centerObject(aoi, 10);

// Display AOI boundary
var empty = ee.Image().byte();
var aoiOutline = empty.paint({
  featureCollection: ee.FeatureCollection(aoi),
  color: 1,
  width: 3
});

// Render on map
Map.addLayer(aoi, {}, 'AOI');




// ==================== SERVER-SIDE DATA COLLECTION FUNCTIONS ====================
/**
 * Returns Image Collection of water mask over the entire year(monsoon season)
 * uses a threshold to mask out water and non water
 * if there are more than 1 images in 2-weeks then it takes the mean of them
 * if there is no image for the biweek(rare) then a complete 0 image is taken(TODO - we can take the prev available image or mean of the prev and next image)
 * Param - year
 **/

function biWeeklyMasks(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, CONFIG.monsoonStart, 1);
  var end = ee.Date.fromYMD(year, CONFIG.monsoonEnd, 30);
  
  var nBiWks = end.difference(start, 'week').divide(2).floor();

  var starts = ee.List.sequence(0, nBiWks.subtract(1))
    .map(function(i) {
      return start.advance(ee.Number(i).multiply(2), 'week');
    });

  var biWeeklyMasksCollection = ee.ImageCollection(starts.map(function(w0) {
    w0 = ee.Date(w0);
    var w1 = w0.advance(2, 'week');

    var col = ee.ImageCollection('COPERNICUS/S1_GRD')
      .filterBounds(aoi)
      .filterDate(w0, w1)
      .filter(ee.Filter.eq('instrumentMode','IW'))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
      .select('VV');

    var isCollectionEmpty = col.size().eq(0);
    var biWeekly = ee.Image(ee.Algorithms.If(isCollectionEmpty, ee.Image().select(), col.mean()));
    var mask = biWeekly.lt(CONFIG.threshold).unmask(0);
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

var yearsEE = ee.List(CONFIG.years);

// Collect bi-weekly masks for all years
var biWeeklyMasksByYear = yearsEE.map(biWeeklyMasks);
print('Data collection initialized for years:', CONFIG.years);




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
function classifyPerennialAndNonWater() {
  var totalBiWeeks = yearsEE.map(function(year) {
    var start = ee.Date.fromYMD(year, CONFIG.monsoonStart, 1);
    var end = ee.Date.fromYMD(year, CONFIG.monsoonEnd, 30);
    return end.difference(start, 'week').divide(2).floor();
  }).reduce(ee.Reducer.sum());

  var waterPresenceSum = ee.Image(0);
  for (var i = 0; i < CONFIG.years.length; i++) {
    var yearCollection = ee.ImageCollection(biWeeklyMasksByYear.get(i));
    waterPresenceSum = waterPresenceSum.add(
      ee.Algorithms.If(
        yearCollection.size().gt(0),
        yearCollection.reduce(ee.Reducer.sum()).unmask(),
        ee.Image(0).toByte()
      )
    );
  }

  var waterPresenceProportion = waterPresenceSum.divide(ee.Image.constant(totalBiWeeks));
  
  var perennialWater = waterPresenceProportion.gte(CONFIG.perennialThreshold).rename('perennial_water');
  var nonWater = waterPresenceProportion.eq(0).rename('non_water');

  var classification = ee.Image(0)
    .where(perennialWater, 1)
    .where(nonWater, 2);

  return classification.clip(aoi);
}

// Initial classification
var perennialAndNonWaterClassification = classifyPerennialAndNonWater();
//visualizing the init classification
// Visualization parameters
var classVis = {
  min: 0,
  max: 2,
  palette: [
    'black', // 0 = neither perennial water nor non-water
    'blue',  // 1 = perennial water
    '#00FF00' // 2 = non-water
  ]
};

// Add to the map
Map.addLayer(perennialAndNonWaterClassification, classVis, 'Perennial & Non-water Classification');
print('Initial classification complete');



/**
*STEP 2...
* creating a function for classifying the unclassified pixels
* This function will be triggered when the inputs are provided and generates a fully classified flood map for 
* the provided period
* 
**/
function processSelectedMask(selectedYear, selectedBiWeek, exportResults) {
  exportResults = exportResults || false;
  
  var selectedYearIndex = CONFIG.years.indexOf(selectedYear);
  var selectedYearCollection = ee.ImageCollection(biWeeklyMasksByYear.get(selectedYearIndex));

  var selectedMask = selectedYearCollection
    .filterMetadata('biweek', 'equals', selectedBiWeek)
    .first();
    
  if (selectedMask) {
    //logging the number of missing images
    // var numExpectedBiWeeks = 9; // Adjust if your season length changes
    // var actualCount = selectedYearCollection.size();
    // var missingCount = ee.Number(numExpectedBiWeeks).subtract(actualCount);
    
    // print('Year:', selectedYear, 
    //       'Expected biweeks:', numExpectedBiWeeks, 
    //       'Available:', actualCount, 
    //       'Missing:', missingCount);
    

    // selectedMask = ee.Image(selectedMask);
    // var selectedBiWeekStartDate = ee.Date.fromYMD(selectedYear, CONFIG.monsoonStart, 1)
    //   .advance(ee.Number(selectedBiWeek).multiply(2), 'week');
    // print('Processing data for bi-week starting:', selectedBiWeekStartDate.format('YYYY-MM-dd'));
      // ----- Added: Compute totalMonsoonBiWeeks same way as old code -----
    var monsoonBiWeeks = ee.List(CONFIG.years).map(function(year) {
      year = ee.Number(year); // ensure it's an ee.Number
      var start = ee.Date.fromYMD(year, CONFIG.monsoonStart, 1);
      var end = ee.Date.fromYMD(year, CONFIG.monsoonEnd, 30);
      return end.difference(start, 'week').divide(2).floor();
    });
    
    var totalMonsoonBiWeeks = ee.Number(monsoonBiWeeks.reduce(ee.Reducer.sum()));

    // -------------------------------------------------------------------

    var waterFrequencyThisBiWeek = ee.ImageCollection(yearsEE.map(function(year) {
      var yearIndex = CONFIG.years.indexOf(year);
      var yearCollection = ee.ImageCollection(biWeeklyMasksByYear.get(yearIndex));
      var biWeekImage = yearCollection.filterMetadata('biweek', 'equals', selectedBiWeek).first();
      return ee.Algorithms.If(biWeekImage, ee.Image(biWeekImage), ee.Image(0).selfMask());
    })).reduce(ee.Reducer.sum()).divide(CONFIG.years.length);


    // ----- Changed: Use totalMonsoonBiWeeks instead of selectedYearCollection.size() -----
    var waterFrequencyThisYear = selectedYearCollection.reduce(ee.Reducer.sum())
      .divide(totalMonsoonBiWeeks);
    // -------------------------------------------------------------------------------------

    var weeklyNonWaterMask = selectedMask.eq(0).rename('weekly_non_water');
    var unclassifiedPixels = perennialAndNonWaterClassification.eq(0);
    var newNonWaterPixels = unclassifiedPixels.and(weeklyNonWaterMask);

    var alreadyClassified = perennialAndNonWaterClassification;
    alreadyClassified = alreadyClassified.where(newNonWaterPixels, 2);

    var seasonalWater = waterFrequencyThisBiWeek.gte(CONFIG.weekFreq)
      .and(waterFrequencyThisYear.lt(CONFIG.yearFreq))
      .rename('seasonal_water')
      .and(alreadyClassified.eq(0));
      
    var floodWater = waterFrequencyThisBiWeek.lt(CONFIG.weekFreq)
      .and(waterFrequencyThisYear.lt(CONFIG.yearFreq))
      .rename('flood_water')
      .and(alreadyClassified.eq(0).and(seasonalWater.not()));
      
    var newPerennialWater = waterFrequencyThisYear.gte(CONFIG.yearFreq)
      .rename('new_seasonal_water')
      .and(alreadyClassified.eq(0).and(seasonalWater.not().and(floodWater.not())));

    var finalClassification = alreadyClassified
      .where(seasonalWater, 3)
      .where(floodWater, 4)
      .where(newPerennialWater, 5);
    
    // Create flood water vectors
    var floodResults = createFloodWaterVectors(finalClassification, selectedYear, selectedBiWeek, exportResults);
    
    // Update map layer
    if (!CONFIG.batchMode) {
      Map.layers().set(1, ui.Map.Layer(finalClassification.clip(aoi), {
        palette: ['000000', '0000FF', '00FF00', 'yellow', 'red', '72A1ED'],
        min: 0,
        max: 5
      }, 'Water Classification'));
    }

    return {
      classification: finalClassification,
      floodVectors: floodResults.vectors,
      floodRaster: floodResults.raster
    };

  } else {
    print('No data found for the selected year and bi-week.');
    if (!CONFIG.batchMode) {
      Map.layers().set(1, ui.Map.Layer(ee.Image().rename('Water Classification'), {}, 'Water Classification', false));
    }
    return null;
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
function createFloodWaterVectors(classificationImage, year, biWeek, exportResults) {
  exportResults = exportResults || false;
  
  if (!CONFIG.batchMode) {
    Map.clear();
    Map.centerObject(aoi, 10);
  }

  var floodMask = classificationImage.eq(3).or(classificationImage.eq(4)).selfMask();
  
  var smudgedMask = floodMask.focal_max({
    radius: 30,
    units: 'meters'
  }).selfMask();

  var floodVectors = smudgedMask.reduceToVectors({
    geometry: aoi,
    scale: 30,
    geometryType: 'polygon',
    labelProperty: 'class',
    maxPixels: 1e13
  });
  
  print('Number of polygons BEFORE area filtering:', floodVectors.size());

  floodVectors = floodVectors.map(function(feature) {
    var area = feature.geometry().area({maxError: 1});
    return feature.set({'area_m2': area});
  });

  var filteredVectors = floodVectors.filter(ee.Filter.gte('area_m2', CONFIG.minAreaSqm));
  
  print('Number of polygons AFTER area filtering:', filteredVectors.size());
  
  // if (!CONFIG.batchMode) {
  //   Map.addLayer(filteredVectors, {color: 'red'}, 'Filtered Flood Polygons');
  // }

  var emptyImage = classificationImage.multiply(0).byte();
  var finalFloodRaster = emptyImage.paint({
    featureCollection: filteredVectors,
    color: 1
  }).clip(aoi);

  if (!CONFIG.batchMode) {
    Map.addLayer(finalFloodRaster.selfMask(), {palette: ['#0000FF']}, 'Final Flood Raster (Masked)');
    Map.addLayer(finalFloodRaster, {min: 0, max: 1, palette: ['black', 'blue']}, 'Final Flood Raster (No Mask)');
  }

  // Export if requested
  if (exportResults || CONFIG.batchMode) {
    exportFloodResults(classificationImage, finalFloodRaster, filteredVectors, year, biWeek);
  }

  return {
    vectors: filteredVectors,
    raster: finalFloodRaster
  };
}


//helper function for exporting images
// ==================== SERVER-SIDE EXPORT FUNCTIONS ====================
function exportFloodResults(classification, floodRaster, floodVectors, year, biWeek) {
  var dateString = year + '_biweek_' + biWeek;
  var regionString = CONFIG.state + '_' + CONFIG.districts.join('_');
  
  if (CONFIG.exportToDrive) {
  
    // Export flood raster
    Export.image.toDrive({
      image: floodRaster,
      description: 'flood_raster_' + regionString + '_' + dateString,
      folder: CONFIG.driveFolder,
      fileNamePrefix: 'flood_raster_' + regionString + '_' + dateString,
      region: aoi,
      scale: classification.projection().nominalScale(),
      crs: classification.projection(),
      maxPixels: 1e13
    });

    // Export flood vectors
    // Export.table.toDrive({
    //   collection: floodVectors,
    //   description: 'flood_vectors_' + regionString + '_' + dateString,
    //   folder: CONFIG.driveFolder,
    //   fileNamePrefix: 'flood_vectors_' + regionString + '_' + dateString,
    //   fileFormat: 'SHP'
    // });
  }

  if (CONFIG.exportToAsset) {
    Export.image.toAsset({
      image: classification,
      description: 'flood_asset_' + regionString + '_' + dateString,
      assetId: CONFIG.assetPath + 'flood_classification_' + regionString + '_' + dateString,
      region: aoi,
      scale: 30,
      maxPixels: 1e13
    });
  }
  
  print('Export tasks submitted for', year, 'bi-week', biWeek);
}

// ==================== BATCH PROCESSING FUNCTIONS ====================
function runBatchProcessing() {
  print('Starting batch processing...');
  CONFIG.batchMode = true;
  
  CONFIG.batchYears.forEach(function(year) {
    CONFIG.batchBiWeeks.forEach(function(biWeek) {
      print('Processing year:', year, 'bi-week:', biWeek);
      processSelectedMask(year, biWeek, true);
    });
  });
  
  print('Batch processing complete. Check Tasks tab for export status.');
  CONFIG.batchMode = false;
}


//***************************UI Elements ****************************************//
// Create UI only if not in batch mode
if (!CONFIG.batchMode) {
  // Main control panel
  var panel = ui.Panel({ style: { width: '350px' } });
  
  // Year selector
  var yearSelect = ui.Select({
    items: CONFIG.years.map(function(year) {
      return { label: String(year), value: year };
    }),
    placeholder: 'Select Year'
  });

  // Bi-week selector  
  var biWeekSelect = ui.Select({
    items: [0,1,2,3,4,5,6,7,8].map(function(week) {
      return { label: 'Bi-Week ' + String(week), value: week };
    }),
    placeholder: 'Select Bi-Week'
  });

  // Interactive display button
  var displayButton = ui.Button({ 
    label: 'Show Bi-Weekly Water',
    onClick: function() {
      var selectedYear = yearSelect.getValue();
      var selectedBiWeek = biWeekSelect.getValue();
      if (selectedYear && selectedBiWeek !== null) {
        processSelectedMask(selectedYear, selectedBiWeek, false);
      } else {
        print('Please select both year and bi-week');
      }
    }
  });

  // Export button
  var exportButton = ui.Button({
    label: 'Export Current Selection',
    onClick: function() {
      var selectedYear = yearSelect.getValue();
      var selectedBiWeek = biWeekSelect.getValue();
      if (selectedYear && selectedBiWeek !== null) {
        processSelectedMask(selectedYear, selectedBiWeek, true);
      } else {
        print('Please select both year and bi-week');
      }
    }
  });

  // Batch processing button
  var batchButton = ui.Button({
    label: 'Run Batch Processing',
    onClick: runBatchProcessing
  });


  // Build panel
  panel.add(ui.Label('Flood Mapping Tool', { fontWeight: 'bold', fontSize: '16px' }));
  panel.add(ui.Label('Select Year and Bi-Week to Display'));
  panel.add(yearSelect);
  panel.add(biWeekSelect);
  panel.add(displayButton);
  panel.add(exportButton);
  panel.add(ui.Label('Batch Processing', { fontWeight: 'bold', fontSize: '14px' }));
  panel.add(ui.Label('Process multiple years/bi-weeks automatically', { fontSize: '12px' }));
  panel.add(batchButton);

  // Legend
  var legendPanel = ui.Panel({
    style: { position: 'bottom-left', padding: '8px 15px' }
  });
  var legendTitle = ui.Label({
    value: 'Water Classification Legend',
    style: { fontWeight: 'bold', fontSize: '16px', margin: '0 0 4px 0' }
  });
  legendPanel.add(legendTitle);
  
  var palette = ['0000FF', '00BB00', 'FFFF00', 'FF0000'];
  var names = ['Perennial Water', 'Non-Water', 'Seasonal Water', 'Flood Water'];
  for (var i = 0; i < palette.length; i++) {
    var colorBox = ui.Label({
      style: { backgroundColor: '#' + palette[i], padding: '8px', margin: '0 0 4px 0' }
    });
    var description = ui.Label({ 
      value: names[i], 
      style: { margin: '0 0 4px 6px' } 
    });
    legendPanel.add(ui.Panel({
      widgets: [colorBox, description],
      layout: ui.Panel.Layout.Flow('horizontal')
    }));
  }
  
  panel.add(legendPanel);
  ui.root.add(panel);
  
  print('UI initialized');
}






