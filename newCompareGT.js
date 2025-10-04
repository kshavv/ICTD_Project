

// Modified version with ROC curve analysis
var CONFIG = {
  // Time parameters
  years: [2018, 2019, 2020, 2021, 2022, 2023, 2024],
  monsoonStart: 5,
  monsoonEnd: 10,
  
  // Classification threshold6
  threshold: -16,
  perennialThreshold: 0.90,
  weekFreq: 0.3,
  yearFreq: 0.3,
  
  // Processing parameters
  minAreaSqm: 100000,
  regionBoundsSize: 0,
  
  // Export parameters
  exportToAsset: false,
  exportToDrive: true,
  assetPath: "users/your_username/flood_maps/",
  driveFolder: "keshav_sparsh/FloodMaps",
  
  // Batch processing mode
  batchMode: true,
  // batchWeekFreq: [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],
  // batchYearFreq: [0.5,0.6,0.7,0.8,0.9,1.0]
  batchWeekFreq: [0.2,0.9],
  batchYearFreq: [0.2]
};

//=========================================================
// Global variables for ROC analysis
var rocResults = [];
var processedCombinations = 0;
var totalCombinations = CONFIG.batchWeekFreq.length * CONFIG.batchYearFreq.length;
var parameterCombinations = [];

// Pre-generate all parameter combinations
function generateParameterCombinations() {
  parameterCombinations = [];
  for (var i = 0; i < CONFIG.batchWeekFreq.length; i++) {
    for (var j = 0; j < CONFIG.batchYearFreq.length; j++) {
      parameterCombinations.push({
        weekFreq: CONFIG.batchWeekFreq[i],
        yearFreq: CONFIG.batchYearFreq[j],
        index: parameterCombinations.length
      });
    }
  }
}

//=========================================================


//========================PROCESSING Ground Truth=================================
//================================================================================
//defining and processing the ground truth
var processedGTImage;
function processGroundTruth(gtImage, minAreaSqm) {

  Map.clear();
  var processingGeom = gtImage.geometry();

  Map.centerObject(processingGeom, 12);
  
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
  


  var rasterFileName = 'GT_flood_raster_';
  var rasterDescription = 'GT_Flood_Raster_';

  return filteredVectors;
}

var gt = gt3;
processedGTImage = processGroundTruth(gt, 300000);

  // Visualize: 0 = black (non-flooded), 1 = green (flooded)
Map.addLayer(processedGTImage, 
              {min: 0, max: 1, palette: ['#000000', '#00FF00']}, 
              'Final GT Raster');
//================================================================================
//================================================================================
//================================================================================


//no prob uptill here ------------------------------------------------------------------- 


//================================================================================
// ==================== SETUP AOI====================
//================================================================================

var aoi = gt.geometry();
Map.centerObject(aoi,10);
Map.addLayer(aoi, {}, 'AOI');

//================================================================================
//================================================================================
//================================================================================



//================================================================================
//===========================DATA COLLECTION=====================================
//================================================================================
/**
* Generates a list of bi-weekly start dates for a given year's monsoon season.
**/
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
**/
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


//================================================================================
//================================================================================
//================================================================================



//================================================================================
/**
**********DATA IS COLLECTED AT THIS POINT. STARTING THE CLASSIFICATION************
**/
//================================================================================

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
function processSelectedMask(selectedYear, selectedBiWeek, yearFreq, weekFreq) {
  var exportResults = false;
  
  var selectedYearMasks = createBiWeeklyMasksForYear(selectedYear);
  var currentWaterMask = ee.Image(selectedYearMasks
    .filter(ee.Filter.eq('biweek_index', selectedBiWeek))
    .first());
  
  
  var seasonStartDate = ee.Date.fromYMD(ee.Number(selectedYear), CONFIG.monsoonStart, 1);
  // Calculate the start date of the specific bi-week.
  var biWeekStartDate = seasonStartDate.advance(ee.Number(selectedBiWeek).multiply(2), 'week');
  
  // Print the date to the console. GEE will automatically format it.
  print('Processing classification for date starting:', biWeekStartDate);
    
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
            (waterFreqThisBiWeek.gte(weekFreq).and(waterFreqThisYear.lt(yearFreq)))
            .or(waterFreqThisBiWeek.lte(weekFreq).and(waterFreqThisYear.gt(yearFreq)))
        );
        
      var floodCondition = unclassified
        .and(currentWaterMask.eq(1))
        .and(waterFreqThisBiWeek.lt(weekFreq));
        
      var newPerennialCondition = unclassified
        .and(currentWaterMask.eq(1))
        .and(waterFreqThisYear.gte(yearFreq)).and(waterFreqThisBiWeek.gte(weekFreq));
  
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

  // Directly create vectors here (server-side)
  // Create vectors server-side, with an empty FeatureCollection fallback
  var floodVectors = ee.FeatureCollection(
    ee.Algorithms.If(
      classificationResult,
      createFloodWaterVectors(finalClassification), // returns ee.FeatureCollection
      ee.FeatureCollection([]) // fallback if no classification (avoid null)
    )
  );

  // Optional visualization if not in batch mode
  if (!CONFIG.batchMode) {
    Map.layers().set(1, ui.Map.Layer(finalClassification.clip(aoi), {
      palette: ['000000', '0000FF', '00FF00', 'FFFF00', 'FF0000'],
      min: 0, max: 4
    }, 'Water Classification'));
  }

  // Return the vector result
  return floodVectors;
}


/**
* Converts flood and seasonal water classes to vector polygons, filters by area, and exports.
*/
function createFloodWaterVectors(classificationImage, year, biWeek, exportResults) {
  exportResults = exportResults || false;
  

  var floodMask = classificationImage.eq(3).or(classificationImage.eq(4)).selfMask();
  // var floodMask = classificationImage.eq(3).selfMask();
  
  
  var smudgedMask = floodMask.focal_max({ radius: 36, units: 'meters' }).selfMask();

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
  // print('Filtered flood/seasonal polygons:', filteredVectors.size());
  
  // // Convert the final filtered vectors back to a raster image
  // var finalFloodRaster = ee.Image(0).byte().paint({
  //   featureCollection: filteredVectors,
  //   color: 1
  // }).clip(aoi.geometry());

  if (!CONFIG.batchMode) {
    Map.addLayer(filteredVectors, {palette: ['#FF5733']}, 'Final Flood Raster');
  }

  // Export if requested
  // if (exportResults || CONFIG.batchMode) {
  //   exportFloodResults(classificationImage, finalFloodRaster, year, biWeek);
  // }

  
    return filteredVectors;
  
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

//================================================================================
//================================================================================
//================================================================================


//================================================================================
//================================================================================
//================================================================================
print('Building historical water masks for all years. This may take a moment...');
var allMasks = buildAllYearsMasks();
var baseInfo = classifyPerennialAndNonWater(allMasks);
var baseClassification = baseInfo.base;
//================================================================================
//================================================================================
//================================================================================




//================================================================================
//====================ADDITIONAL FUNCTIONS FOR ROC CURVE ANALYSIS=================
//================================================================================

function calculateTPRandFPRFromVectors(predictedVectors, gtVectors,weekFreq,yearFreq,callback) {
  // Convert both vector collections to 1/0 rasters
  var predictedRaster = ee.Image(0).byte().paint({
    featureCollection: predictedVectors,
    color: 1
  }).clip(aoi);

  var gtRaster = ee.Image(0).byte().paint({
    featureCollection: gtVectors,
    color: 1
  }).clip(aoi);

  // Compute TP, FP, FN, TN masks
  var truePositives  = predictedRaster.and(gtRaster);
  var falsePositives = predictedRaster.and(gtRaster.not());
  var falseNegatives = gtRaster.and(predictedRaster.not());
  var trueNegatives  = predictedRaster.not().and(gtRaster.not());

  // Combine all in one reduceRegion call
  var combined = ee.Image([
    truePositives.rename('TP'),
    falsePositives.rename('FP'),
    falseNegatives.rename('FN'),
    trueNegatives.rename('TN')
  ]);

  combined.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e13
  }).evaluate(function(result, error) {
    if (error) {
      print('Error calculating metrics for weekFreq:', weekFreq, 'yearFreq:', yearFreq, error);
      callback(null);
    } else {
      var tp = result.TP || 0;
      var fp = result.FP || 0;
      var fn = result.FN || 0;
      var tn = result.TN || 0;

      // Avoid division by zero
      var tpr = (tp + fn) > 0 ? tp / (tp + fn) : 0;
      var fpr = (fp + tn) > 0 ? fp / (fp + tn) : 0;

      callback({
        TPR: tpr,
        FPR: fpr,
        weekFreq: weekFreq,
        yearFreq: yearFreq,
        TP: tp,
        FP: fp,
        FN: fn,
        TN: tn
      });
    }
  });
}
// Sequential processing function to avoid race conditions
function processNextCombination(combinationIndex) {
  if (combinationIndex >= parameterCombinations.length) {
    // All combinations processed, generate ROC curve
    print('All combinations processed. Generating ROC curve...');
    generateROCCurve();
    return;
  }
  
  var combo = parameterCombinations[combinationIndex];
  var weekFreq = combo.weekFreq;
  var yearFreq = combo.yearFreq;
  
  print('Processing combination:', combinationIndex + 1, 'of', totalCombinations, 
        '- WeekFreq:', weekFreq, 'YearFreq:', yearFreq);
  
  
  //processSelectedMasks(selectedYear,selectedWeek,yearFreq,weekFreq)
  var result = processSelectedMask(2018, 3, yearFreq, weekFreq);
  var layerName = 'Flood polygons' + '_Yfreq' + yearFreq + '_Wfreq' + weekFreq;

  Map.addLayer(result, {}, layerName);
  
  if (result) {
      
       calculateTPRandFPRFromVectors(result, processedGTImage,weekFreq,yearFreq, function(metrics) {
      if (metrics) {
        rocResults.push(metrics);
        print('Combination', combinationIndex + 1, 'complete - TPR:', metrics.TPR.toFixed(4), 'FPR:', metrics.FPR.toFixed(4));
      } else {
        print('Failed to calculate metrics for combination', combinationIndex + 1);
      }
      
      // Process next combination only after current one is complete
      processNextCombination(combinationIndex + 1);
    });
    
  } else {
    print('No results for combination', combinationIndex + 1, '- WeekFreq:', weekFreq, 'YearFreq:', yearFreq);
    // Continue to next combination immediately if no processing needed
    processNextCombination(combinationIndex + 1);
  }
  
  if(result && result.floodResults){//update thsi once we are getting results
      // processNextCombination(combinationIndex + 1);  
  }
}

// Batch processing function
function runBatchProcessing() {
  print('Starting ROC analysis...');
  print('Total combinations to process:', totalCombinations);
  
  // Generate all parameter combinations first
  generateParameterCombinations();
  
  print('Parameter combinations generated. Starting sequential processing...');
  
  // Start processing from the first combination
  processNextCombination(0);
}

// Function to generate and display ROC curve
function generateROCCurve() {
  print('Generating ROC curve with', rocResults.length, 'data points...');
  
  // Sort results by FPR for proper curve plotting
  rocResults.sort(function(a, b) {
    return a.FPR - b.FPR;
  });
  
  // Print results in CSV format for easy copying
  print('ROC Results (CSV format):');
  print('WeekFreq,YearFreq,TPR,FPR,TP,FP,FN,TN');
  
  for (var i = 0; i < rocResults.length; i++) {
    var r = rocResults[i];
    print(r.weekFreq + ',' + r.yearFreq + ',' + r.TPR.toFixed(6) + ',' + 
          r.FPR.toFixed(6) + ',' + r.TP + ',' + r.FP + ',' + r.FN + ',' + r.TN);
  }
  
  // Calculate AUC using trapezoidal rule
  var auc = 0;
  for (var i = 1; i < rocResults.length; i++) {
    var dx = rocResults[i].FPR - rocResults[i-1].FPR;
    var avgY = (rocResults[i].TPR + rocResults[i-1].TPR) / 2;
    auc += dx * avgY;
  }
  
  print('Calculated AUC:', auc.toFixed(4));
  
  // Create a simple chart using GEE's built-in charting
  var chartData = [];
  for (var i = 0; i < rocResults.length; i++) {
    chartData.push([rocResults[i].FPR, rocResults[i].TPR]);
  }
  
  // Add diagonal line points for reference
  chartData.push([0, 0]);
  chartData.push([1, 1]);
  
  print('ROC Curve Data Points (FPR, TPR):');
  for (var i = 0; i < chartData.length - 2; i++) {
    print('(' + chartData[i][0].toFixed(4) + ', ' + chartData[i][1].toFixed(4) + ')');
  }
  
  print('Analysis complete! Copy the CSV data above to create ROC curve visualization.');
  print('Best performing combination (highest TPR-FPR):');
  
  var bestMetric = -Infinity;
  var bestCombination = null;
  for (var i = 0; i < rocResults.length; i++) {
    var metric = rocResults[i].TPR - rocResults[i].FPR; // Youden's index
    if (metric > bestMetric) {
      bestMetric = metric;
      bestCombination = rocResults[i];
    }
  }
  
  if (bestCombination) {
    print('Best: WeekFreq =', bestCombination.weekFreq, 
          ', YearFreq =', bestCombination.yearFreq,
          ', TPR =', bestCombination.TPR.toFixed(4),
          ', FPR =', bestCombination.FPR.toFixed(4),
          ', Youden Index =', bestMetric.toFixed(4));
  }
}

// Start the sequential processing
print('Initializing ROC analysis with', totalCombinations, 'parameter combinations...');
print('This will process combinations sequentially to avoid race conditions.');
runBatchProcessing();
