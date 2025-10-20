

// ==================== CONFIGURATION PARAMETERS ====================
var CONFIG = {
  // Time parameters
  years: [2018, 2019, 2020, 2021, 2022, 2023, 2024],
  monsoonStart: 5, // Inclusive start month
  monsoonEnd: 10,   // Inclusive end month
  
  // Classification thresholds
  threshold: -16,           
  perennialThreshold: 0.9,  
  weekFreq: 0.6,            
  yearFreq: 0.9,            
  
  // Processing parameters
  minAreaSqm: 100000,
  regionBoundsSize: 0, // Set to 0 to use full AOI bounds
  
  // Export parameters
  exportToAsset: false,
  exportToDrive: true,
  // assetPath: "users/your_username/flood_maps/",
  driveFolder: "GEE_Flood_Exports",
  
  // Batch processing mode
  batchMode: false, // Set to true for automated batch processing
  batchYears: [2023, 2024], // Years to process in batch mode
  batchBiWeeks: [3, 4, 5] // Bi-weeks to process in batch mode
};


//====================PROCESSING THE GT IMAGE ======================
var processedGTImage;
function processGroundTruth(gtImage, minAreaSqm) {

  Map.clear();
  var processingGeom = gtImage.geometry();

  // Map.centerObject(processingGeom, 12);
  
  // Set default export scale if not provided
  var exportScale = 30;

  // Isolate flood pixels (value is 1 for ground truth).
  var floodMask = gtImage.eq(1).selfMask();
  
  // The 'focal_max' (smudging) step is skipped for ground truth data
  // to preserve its original accuracy.
  // Step 2: Vectorize the mask to create polygons from the raster data.
  var floodVectors = floodMask.reduceToVectors({
    geometry: processingGeom, // Use the image's geometry
    scale: exportScale, // Use the export scale for vectorization
    geometryType: 'polygon',
    labelProperty: 'class',
    maxPixels: 1e13
  });
  
  // Step 3: Compute the area of each polygon and store it as a property.
  floodVectors = floodVectors.map(function(feature) {
    var area = feature.geometry().area({maxError: 1});
    return feature.set({'area_m2': area});
  });

  // Step 4: Filter out polygons that are smaller than the specified minimum area.
  var filteredVectors = floodVectors.filter(ee.Filter.gte('area_m2', minAreaSqm));

  
  // Step 5: Convert the filtered vectors back into a raster image.
  // Rasterize filteredVectors: flood areas get value 1, background is 0
  var finalGTRaster = ee.Image(0).byte().paint({
    featureCollection: filteredVectors,
    color: 1
  }).clip(processingGeom);
  
  // Visualize: 0 = black (non-flooded), 1 = green (flooded)
  Map.addLayer(filteredVectors, 
              {min: 0, max: 1, palette: ['#000000', '#00FF00']}, 
              'Final GT Raster');

  var rasterFileName = 'GT_flood_raster_';
  var rasterDescription = 'GT_Flood_Raster_';

  return finalGTRaster;
}
processedGTImage = processGroundTruth(gt, 300000);


// ==================== REGION SETUP ====================

  var districts = ee.FeatureCollection("FAO/GAUL/2015/level2");
  var myPredictionRegion = districts.filter(ee.Filter.and(
    ee.Filter.eq("ADM1_NAME", "Kerala"),
      ee.Filter.or(
      // ee.Filter.eq("ADM2_NAME", "Ernakulam"),
      ee.Filter.eq("ADM2_NAME", "Kottayam"),
      ee.Filter.eq("ADM2_NAME", "Alappuzha"),
      ee.Filter.eq("ADM2_NAME", "Pattanamtitta")
    )
  ));
  
  Map.addLayer(myPredictionRegion,{},'my_predicred_region');


function setupRegion() {


  var aoi;
  if (CONFIG.regionBoundsSize) { //clip the boundary based on the clipping box
    aoi = region.geometry().centroid().buffer(CONFIG.regionBoundsSize / 2).bounds();
    print('Using centered export box of size (m):', CONFIG.regionBoundsSize);
  } else {
    aoi = myPredictionRegion; //return the entire region
    // aoi = processedGTImage.geometry(); 
    // aoi = gt;
    print('Using default AOI bounds for export.');
  }
  
  return aoi;
}

// Initialize AOI
var aoi = setupRegion();
Map.centerObject(aoi);
Map.addLayer(aoi, {}, 'AOI');


// ==================== DATA COLLECTION (LOGIC FROM SCRIPT 2) ====================

/**
 * Generates a list of bi-weekly start dates for a given year's monsoon season.
 */
function listBiWeekStarts(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, CONFIG.monsoonStart, 1);
  var end = ee.Date.fromYMD(year, CONFIG.monsoonEnd, 1).advance(1, 'month').advance(-1, 'day');
  var nBi = end.difference(start, 'week').divide(2).ceil();
  return ee.List.sequence(0, nBi.subtract(1)).map(function(i) {
    return start.advance(ee.Number(i).multiply(2), 'week');
  });
}

/**
 * Retrieves Sentinel-1 imagery for a specific period and clips it to the AOI.
 */
function s1CollectionForPeriod(start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .select(['VV'])
    // .map(function(img) { return img.clip(aoi.geometry()); });
    .map(function(img) { return img.clip(aoi); });
}

/**
 * Creates bi-weekly water masks for a single year.
 * If a bi-week has no S1 images, it creates a fully masked image (a data gap).
 * @param {ee.Number|number} year The year to process.
 * @return {ee.ImageCollection} The collection of bi-weekly water masks.
 */
function createBiWeeklyMasksForYear(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, CONFIG.monsoonStart, 1);
  var starts = listBiWeekStarts(year);

  return ee.ImageCollection(
    starts.map(function(w0) {
      w0 = ee.Date(w0);
      var w1 = w0.advance(2, 'week');

      var s1col = s1CollectionForPeriod(w0, w1);
      var hasData = s1col.size().gt(0);

      // If data exists, compute the mean. Otherwise, create a masked image.
      var biweeklyMean = ee.Algorithms.If(
        hasData, 
        s1col.mean(), 
        ee.Image(0).updateMask(ee.Image(0)) // Properly masked image for gaps
      );
      biweeklyMean = ee.Image(biweeklyMean);

      var mask = biweeklyMean.lt(CONFIG.threshold).rename('water');
      var biWkNum = w0.difference(start, 'week').divide(2).floor();

      return mask.set({
        'year': year,
        'biweek_index': biWkNum,
        'system:time_start': w0.millis()
      }).clip(aoi); //debug (original .clip(aoi.geometry())) -- key::geometry
    })
  );
}

/**
 * Builds a single ImageCollection containing all bi-weekly masks from all years.
 */
function buildAllYearsMasks() {
  var allImagesList = ee.List(CONFIG.years).map(function(y) {
    // For each year, get the collection and convert it to a list of images.
    // 100 is a safe upper limit for the number of bi-weeks in a monsoon season.
    return createBiWeeklyMasksForYear(y).toList(100);
  }).flatten(); // Flatten the list of lists into a single list of images.
  
  // Create a single collection from the flat list of all images.
  return ee.ImageCollection(allImagesList);
}


// ==================== CLASSIFICATION (LOGIC FROM SCRIPT 2) ====================

/**
 * STEP 1: Classifies pixels as perennial water or permanent non-water
 * based on the entire historical dataset.
 * @param {ee.ImageCollection} allMasksAllYears - Collection from buildAllYearsMasks().
 * @return {Object} An object containing the base classification and frequency info.
 */
function classifyPerennialAndNonWater(allMasksAllYears) {
  // Count how many valid (unmasked) observations exist for each pixel
  var validCount = allMasksAllYears
    .map(function(img) { return img.mask().rename('m'); })
    .reduce(ee.Reducer.sum())
    .unmask(0); // Ensure a valid image even if there are no observations

  // Sum the water presence (where water=1)
  var waterPresenceSum = allMasksAllYears.reduce(ee.Reducer.sum())
    .unmask(0);

  // Calculate frequency: sum of water / number of valid observations
  var freq = waterPresenceSum.divide(validCount.add(1e-10)) // avoid division by zero
    .updateMask(validCount.gt(0));

  var perennial = freq.gte(CONFIG.perennialThreshold);
  var nonWater = freq.eq(0);

  // Create base classification: 1 for perennial, 2 for non-water, 0 for others
  var baseCls = ee.Image(0)
    .where(perennial, 1)
    .where(nonWater, 2)
    .rename('classification');
    
  return {
    base: baseCls.clip(aoi), //debug (original .clip(aoi.geometry())) -- key::geometry
    freq: freq
  };
}

/**
 * STEP 2: Performs the final classification for a selected year and bi-week.
 * Classifies the remaining pixels into seasonal, flood, or temporary non-water.
 * This is the main function triggered by the UI.
 */
function processSelectedMask(selectedYear, selectedBiWeek, exportResults) {
  exportResults = exportResults || false;
  
  var selectedYearMasks = createBiWeeklyMasksForYear(selectedYear);
  
  // This may be null if no image matches the filter
  var currentWaterMask = ee.Image(selectedYearMasks
    .filter(ee.Filter.eq('biweek_index', selectedBiWeek))
    .first());
    
  // Use ee.Algorithms.If to handle cases where no image is found server-side.
  // This entire block remains on the server.
  var classificationResult = ee.Algorithms.If(
    currentWaterMask, // Condition: checks if currentWaterMask is not null
    
    // --- IF TRUE: An image was found, so run the full classification ---
    (function() {
      var validCountYear = selectedYearMasks.map(function(img) { return img.mask(); })
        .reduce(ee.Reducer.sum()).unmask(0);
      var waterSumYear = selectedYearMasks.reduce(ee.Reducer.sum()).unmask(0);
      var waterFreqThisYear = waterSumYear.divide(validCountYear.add(1e-10))
        .updateMask(validCountYear.gt(0));

      var biweekAcrossYears = allMasks.filter(ee.Filter.eq('biweek_index', selectedBiWeek));
      var validCountBi = biweekAcrossYears.map(function(i) { return i.mask(); })
        .reduce(ee.Reducer.sum()).unmask(0);
      var waterSumBi = biweekAcrossYears.reduce(ee.Reducer.sum()).unmask(0);
      var waterFreqThisBiWeek = waterSumBi.divide(validCountBi.add(1e-10))
        .updateMask(validCountBi.gt(0));
      
      var unclassified = baseClassification.eq(0);
      
      var seasonalCondition = unclassified
        .and(currentWaterMask.eq(1))
        .and
        (
            (waterFreqThisBiWeek.gte(CONFIG.weekFreq).and(waterFreqThisYear.lt(CONFIG.yearFreq)))
            .or(waterFreqThisBiWeek.lte(CONFIG.weekFreq).and(waterFreqThisYear.gt(CONFIG.yearFreq)))
        );
        
      var floodCondition = unclassified
        .and(currentWaterMask.eq(1))
        .and(waterFreqThisBiWeek.lt(CONFIG.weekFreq));
        
      var newPerennialCondition = unclassified
        .and(currentWaterMask.eq(1))
        .and(waterFreqThisYear.gte(CONFIG.yearFreq)).and(waterFreqThisBiWeek.gte(CONFIG.weekFreq));
  
      var temporaryNonWater = unclassified.and(currentWaterMask.eq(0));
  
      return baseClassification
        .where(temporaryNonWater, 2)
        .where(seasonalCondition, 3)
        .where(floodCondition, 4)
        .where(newPerennialCondition, 1);
    })(),
    
    // --- IF FALSE: No image was found, return a placeholder ---
    null // Return null to signify no data
  );

  // Cast the server-side result to an ee.Image. If it was null, this will be a
  // computed object that represents null.
  var finalClassification = ee.Image(classificationResult);

  // Now, evaluate a property of the result (like its band names) to check if it's a
  // valid image on the client side.
  finalClassification.bandNames().evaluate(function(bands, error) {
    if (error) {
      print('An error occurred during classification:', error);
      return;
    }

    // If the 'bands' list exists and is not empty, the classification was successful.
    if (bands && bands.length > 0) {
      var floodResults = createFloodWaterVectors(finalClassification, selectedYear, selectedBiWeek, exportResults);
      if (!CONFIG.batchMode) {
        Map.layers().set(1, ui.Map.Layer(finalClassification.clip(aoi), {//debug (original .clip(aoi.geometry())) -- key::geometry
          palette: ['000000', '0000FF', '00FF0000', 'FFFF00', 'FF0000'],
          min: 0, max: 4
        }, 'Water Classification'));
      }
    } else {
      // If 'bands' is null or empty, it means no image was created (no data).
      print('No data found for the selected year and bi-week.');
      if (!CONFIG.batchMode) {
        Map.layers().set(1, ui.Map.Layer(ee.Image(), {}, 'No Data', false));
      }
    }
  });
}



// ==================== VECTORIZATION & EXPORT (FROM SCRIPT 1) ====================

/**
 * Converts flood and seasonal water classes to vector polygons, filters by area, and exports.
 */
function createFloodWaterVectors(classificationImage, year, biWeek, exportResults) {
  exportResults = exportResults || false;
  

  var floodMask = classificationImage.eq(3).or(classificationImage.eq(4)).selfMask();
  // var floodMask = classificationImage.eq(3).selfMask();
  
  
  var smudgedMask = floodMask.focal_max({ radius: 50, units: 'meters' }).selfMask();

  var floodVectors = smudgedMask.reduceToVectors({
    geometry: aoi, //debug (original .clip(aoi.geometry())) -- key::geometry
    scale: 150,
    geometryType: 'polygon',
    maxPixels: 1e13
  });
  
  // Map over the vectors to calculate the area of each polygon
  floodVectors = floodVectors.map(function(feature) {
    return feature.set({'area_m2': feature.geometry().area({maxError: 1})});
  });

  // Filter the polygons by the minimum area threshold
  var filteredVectors = floodVectors.filter(ee.Filter.gte('area_m2', CONFIG.minAreaSqm));
  print('Filtered flood/seasonal polygons:', filteredVectors.size());
  
  // // Convert the final filtered vectors back to a raster image
  // var finalFloodRaster = ee.Image(0).byte().paint({
  //   featureCollection: filteredVectors,
  //   color: 1
  // }).clip(aoi.geometry());

  if (!CONFIG.batchMode) {
    Map.addLayer(filteredVectors, {palette: ['#FF5733']}, 'Final Flood Raster');
  }

  // Export if requested
  if (exportResults || CONFIG.batchMode) {
    exportFloodResults(classificationImage, finalFloodRaster, year, biWeek);
  }

  return {
    vectors: filteredVectors
  };
}

/**
 * Handles the export tasks to Google Drive or Earth Engine Assets.
 */
function exportFloodResults(classification, floodRaster, year, biWeek) {
  var dateString = year + '_biweek_' + biWeek;
  
  if (CONFIG.exportToDrive) {
    Export.image.toDrive({
      image: floodRaster,
      description: 'flood_raster_' + dateString,
      folder: CONFIG.driveFolder,
      fileNamePrefix: 'flood_raster_' + dateString,
      region: aoi.geometry(),
      scale: 30,
      maxPixels: 1e13
    });
  }

  if (CONFIG.exportToAsset) {
    Export.image.toAsset({
      image: classification,
      description: 'flood_classification_asset_' + dateString,
      assetId: CONFIG.assetPath + 'flood_classification_' + dateString,
      region: aoi.geometry(),
      scale: 30,
      maxPixels: 1e13
    });
  }
  
  print('Export tasks submitted for', year, 'bi-week', biWeek);
}


// ==================== BATCH PROCESSING (FROM SCRIPT 1) ====================
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


// ==================== INITIALIZATION & UI ====================

// --- Main Data Loading ---
print('Building historical water masks for all years. This may take a moment...');
var allMasks = buildAllYearsMasks();
var baseInfo = classifyPerennialAndNonWater(allMasks);
var baseClassification = baseInfo.base;

Map.addLayer(baseClassification, {
  min: 0, max: 2, palette: ['black', 'blue', 'green']
}, 'Base Classification (Perennial/Non-Water)', false);
print('Base classification complete.');


// --- UI Elements ---
if (!CONFIG.batchMode) {
  var panel = ui.Panel({ style: { width: '350px' } });
  
  var yearSelect = ui.Select({
    items: CONFIG.years.map(String).map(function(y) { return {label: y, value: parseInt(y, 10)}; }),
    placeholder: 'Select Year'
  });

  var biWeekSelect = ui.Select({
    items: ee.List.sequence(0, 12).getInfo().map(function(w) { return {label: 'Bi-Week ' + w, value: w}; }),
    placeholder: 'Select Bi-Week'
  });

  var displayButton = ui.Button({
    label: 'Generate Flood Map',
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

  var batchButton = ui.Button({
    label: 'Run Batch Processing',
    onClick: runBatchProcessing
  });

  panel.add(ui.Label('Flood Mapping Tool', { fontWeight: 'bold', fontSize: '16px' }));
  panel.add(ui.Label('Select Year and Bi-Week to Display'));
  panel.add(yearSelect);
  panel.add(biWeekSelect);
  panel.add(displayButton);
  panel.add(exportButton);
  panel.add(ui.Label('Batch Processing', { fontWeight: 'bold', fontSize: '14px', margin: '10px 0 0 0' }));
  panel.add(batchButton);

  var legendPanel = ui.Panel({
    style: { position: 'bottom-left', padding: '8px 15px' }
  });
  legendPanel.add(ui.Label({
    value: 'Legend', style: { fontWeight: 'bold', fontSize: '16px', margin: '0 0 4px 0' }
  }));
  
  var palette = ['0000FF', '00FF00', 'FFFF00', 'FF0000'];
  var names = ['1: Perennial Water', '2: Non-Water', '3: Seasonal Water', '4: Flood Water'];
  for (var i = 0; i < palette.length; i++) {
    var colorBox = ui.Label({ style: { backgroundColor: '#' + palette[i], padding: '8px', margin: '0 0 4px 0' } });
    var description = ui.Label({ value: names[i], style: { margin: '0 0 4px 6px' } });
    legendPanel.add(ui.Panel({
      widgets: [colorBox, description],
      layout: ui.Panel.Layout.Flow('horizontal')
    }));
  }
  
  ui.root.add(panel);
  Map.add(legendPanel);
  print('UI initialized.');
}

