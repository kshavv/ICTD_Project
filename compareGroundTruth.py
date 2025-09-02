

// Modified version with ROC curve analysis
var CONFIG = {
  // Time parameters
  years: [2018, 2019, 2020, 2021, 2022, 2023, 2024],
  monsoonStart: 5,
  monsoonEnd: 10,
  
  // Classification threshold6
  threshold: -16,
  perennialThreshold: 0.95,
  weekFreq: 0.3,
  yearFreq: 0.3,
  
  // Processing parameters
  minAreaSqm: 50000,
  regionBoundsSize: 0,
  
  // Export parameters
  exportToAsset: false,
  exportToDrive: true,
  assetPath: "users/your_username/flood_maps/",
  driveFolder: "keshav_sparsh/FloodMaps",
  
  // Batch processing mode
  batchMode: true,

  batchWeekFreq: [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],
  batchYearFreq: [0.8]
  // batchWeekFreq: [0.2,0.3,0.6,0.7],
  // batchYearFreq: [0.2,0.3,0.6,0.7]
};

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

//defining and processing the ground truth
var processedGTImage;
function processGroundTruth(gtImage, minAreaSqm) {
  // Clear the map of any previous layers before running.
  Map.clear();
  var processingGeom = gtImage.geometry();
  
  // Center the map on the image's geometry.
  Map.centerObject(processingGeom, 10);
  
  // Set default export scale if not provided
  var exportScale = 30;


  // --- MODIFIED STEP 1 ---
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
  
  print('Number of GT polygons after filtering:', filteredVectors.size());
  
  
  // Step 5: Convert the filtered vectors back into a raster image.
  // Rasterize filteredVectors: flood areas get value 1, background is 0
  var finalGTRaster = ee.Image(0).byte().paint({
    featureCollection: filteredVectors,
    color: 1
  }).clip(processingGeom);
  
  // Visualize: 0 = black (non-flooded), 1 = green (flooded)
  Map.addLayer(finalGTRaster, 
              {min: 0, max: 1, palette: ['#000000', '#00FF00']}, 
              'Final GT Raster');

  var rasterFileName = 'GT_flood_raster_';
  var rasterDescription = 'GT_Flood_Raster_';

  Export.image.toDrive({
    image: finalGTRaster,
    description: rasterDescription,
    folder: 'keshav_sparsh/FloodMaps',
    fileNamePrefix: rasterFileName,
    region: processingGeom.bounds(), // Use the bounds of the image's geometry for export
    scale: exportScale,
    crs: gtImage.projection()
  });

  // Step 7: Return the final raster image.
  return finalGTRaster;
}
processedGTImage = processGroundTruth(gt1, 500000);


// ==================== SETUP AOI====================
// Define AOI from processed GT image
var aoi = processedGTImage.geometry(); // Uses exact extent of GT raster
Map.centerObject(aoi, 10);



//************************collect the biweekly images*********************************//
function biWeeklyMasks(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, CONFIG.monsoonStart, 1);
  var end = ee.Date.fromYMD(year, CONFIG.monsoonEnd, 30);
  
  var nBiWks = end.difference(start, 'week').divide(2).floor();

  var starts = ee.List.sequence(0, nBiWks.subtract(1))
    .map(function(i) {
      return start.advance(ee.Number(i).multiply(2), 'week');
    });

  // Initial state for the iteration: an empty collection and an empty "last valid image".
  var initialState = ee.Dictionary({
    'collection': ee.ImageCollection([]),
    'lastImage': ee.Image().select() // An image with no bands is our placeholder for "none yet"
  });

  // The function that will be applied to each bi-week start date.
  var accumulate = function(w0, previousState) {
    // Cast the inputs to the correct types.
    w0 = ee.Date(w0);
    previousState = ee.Dictionary(previousState);
    var collection = ee.ImageCollection(previousState.get('collection'));
    var lastImage = ee.Image(previousState.get('lastImage'));

    var w1 = w0.advance(2, 'week');

    // Get the S1 collection for the current bi-week.
    var col = ee.ImageCollection('COPERNICUS/S1_GRD')
      .filterBounds(aoi)
      .filterDate(w0, w1)
      .filter(ee.Filter.eq('instrumentMode','IW'))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
      .select('VV');

    var isCollectionEmpty = col.size().eq(0);
    
    // If the collection is NOT empty, the current image is the mean.
    // If it IS empty, we carry forward the last valid image.
    var newMean = col.mean();
    var imageToProcess = ee.Image(ee.Algorithms.If(isCollectionEmpty, lastImage, newMean));

    // **Edge Case Handler**: If imageToProcess has no bands, it means this is the
    // very first bi-week and it was empty. In this case, we fall back to a zero image.
    // Otherwise, we use the valid new image or the valid carried-over image.
    var finalImage = ee.Image(ee.Algorithms.If(
        imageToProcess.bandNames().size().eq(0),
        ee.Image(0), // Fallback for the very first step
        imageToProcess
    ));

    // Create the mask from the determined image (new, carried-over, or fallback).
    var mask = finalImage.lt(CONFIG.threshold).unmask(0);
    var biWkNum = w0.difference(start, 'week').divide(2).floor();

    var newMask = mask
      .rename("water")
      .set({
        'year': year,
        'biweek': biWkNum,
        'system:time_start': w0.millis(),
        'has_data': ee.Algorithms.If(isCollectionEmpty, 0, 1) // Note: this now means "had new data this period"
      })
      .clip(aoi);

    // Add the new mask to our collection.
    var newCollection = collection.merge(ee.ImageCollection([newMask]));
    
    // Update the "lastImage" state ONLY if we got new data this period.
    var newLastImage = ee.Image(ee.Algorithms.If(isCollectionEmpty, lastImage, newMean));

    // Return the new state for the next iteration.
    return ee.Dictionary({
      'collection': newCollection,
      'lastImage': newLastImage
    });
  };

  // Run the iteration over all bi-week start dates.
  var finalState = ee.Dictionary(starts.iterate(accumulate, initialState));
  
  // Extract the final collection from the state dictionary.
  return ee.ImageCollection(finalState.get('collection'));
}

var yearsEE = ee.List(CONFIG.years);

// Collect bi-weekly masks for all years
var biWeeklyMasksByYear = yearsEE.map(biWeeklyMasks);
print('Data collection initialized for years:', CONFIG.years);




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


/**
*STEP 2...
* creating a function for classifying the unclassified pixels
* This function will be triggered when the inputs are provided and generates a fully classified flood map for 
* the provided period
* 
**/
// Modified processSelectedMask function to return metrics directly
function processSelectedMask(selectedYear, selectedBiWeek, yearFreq, weekFreq) {
  var selectedYearIndex = CONFIG.years.indexOf(selectedYear);
  var selectedYearCollection = ee.ImageCollection(biWeeklyMasksByYear.get(selectedYearIndex));
  
  var selectedMask = selectedYearCollection
    .filterMetadata('biweek', 'equals', selectedBiWeek)
    .first();
    
  // Calculate the start date of the bi-week being processed
  var start = ee.Date.fromYMD(ee.Number(selectedYear), CONFIG.monsoonStart, 1);
  var biWeekStartDate = start.advance(ee.Number(selectedBiWeek).multiply(2), 'week');
  print('Processing date:', biWeekStartDate);
    
  if (selectedMask) {
    var monsoonBiWeeks = ee.List(CONFIG.years).map(function(year) {
      year = ee.Number(year);
      var start = ee.Date.fromYMD(year, CONFIG.monsoonStart, 1);
      var end = ee.Date.fromYMD(year, CONFIG.monsoonEnd, 30);
      return end.difference(start, 'week').divide(2).floor();
    });
    
    var totalMonsoonBiWeeks = ee.Number(monsoonBiWeeks.reduce(ee.Reducer.sum()));

    var waterFrequencyThisBiWeek = ee.ImageCollection(yearsEE.map(function(year) {
      var yearIndex = CONFIG.years.indexOf(year);
      var yearCollection = ee.ImageCollection(biWeeklyMasksByYear.get(yearIndex));
      var biWeekImage = yearCollection.filterMetadata('biweek', 'equals', selectedBiWeek).first();
      return ee.Algorithms.If(biWeekImage, ee.Image(biWeekImage), ee.Image(0).selfMask());
    })).reduce(ee.Reducer.sum()).divide(CONFIG.years.length);

    var waterFrequencyThisYear = selectedYearCollection.reduce(ee.Reducer.sum())
      .divide(totalMonsoonBiWeeks);

    var weeklyNonWaterMask = selectedMask.eq(0).rename('weekly_non_water');
    var unclassifiedPixels = perennialAndNonWaterClassification.eq(0);
    var newNonWaterPixels = unclassifiedPixels.and(weeklyNonWaterMask);

    var alreadyClassified = perennialAndNonWaterClassification;
    alreadyClassified = alreadyClassified.where(newNonWaterPixels, 2);

    var seasonalWater = waterFrequencyThisBiWeek.gte(weekFreq)
      .and(waterFrequencyThisYear.lt(yearFreq))
      .rename('seasonal_water')
      .and(alreadyClassified.eq(0));
      
    var floodWater = waterFrequencyThisBiWeek.lt(weekFreq)
      .and(waterFrequencyThisYear.lt(yearFreq))
      .rename('flood_water')
      .and(alreadyClassified.eq(0).and(seasonalWater.not()));
      
    var newPerennialWater = waterFrequencyThisYear.gte(yearFreq)
      .rename('new_seasonal_water')
      .and(alreadyClassified.eq(0).and(seasonalWater.not().and(floodWater.not())));

    var finalClassification = alreadyClassified
      .where(seasonalWater, 3)
      .where(floodWater, 4)
      .where(newPerennialWater, 5);
    
    // var floodResults = createFloodWaterVectors(finalClassification);
    var floodResults = finalClassification.eq(4);
    var floodMapName = 'Flood-Week_' + yearFreq +'_'+ weekFreq;
    
    // Map.addLayer(floodResults, {min: 0, max: 1, palette: ['black', 'blue']}, floodMapName);
    Map.layers().set(1, ui.Map.Layer(finalClassification.clip(aoi), {
      palette: ['000000', '0000FF', '00FF00', 'yellow', 'red', '72A1ED'],
      min: 0,
      max: 5
    }, floodMapName));
    
    return {
      floodResults: floodResults,
      weekFreq: weekFreq,
      yearFreq: yearFreq
    };
  } else {
    return null;
  }
}

function createFloodWaterVectors(classificationImage) {
  var floodMask = classificationImage.eq(3).or(classificationImage.eq(4)).selfMask();
  
  var smudgedMask = floodMask.focal_max({
    radius: 25,
    units: 'meters'
  }).selfMask();

  var floodVectors = smudgedMask.reduceToVectors({
    geometry: aoi,
    scale: 30,
    geometryType: 'polygon',
    labelProperty: 'class',
    maxPixels: 1e13
  });

  floodVectors = floodVectors.map(function(feature) {
    var area = feature.geometry().area({maxError: 1});
    return feature.set({'area_m2': area});
  });

  var filteredVectors = floodVectors.filter(ee.Filter.gte('area_m2', CONFIG.minAreaSqm));

  var emptyImage = classificationImage.multiply(0).byte();
  var finalFloodRaster = emptyImage.paint({
    featureCollection: filteredVectors,
    color: 1
  }).clip(aoi);
  

  return finalFloodRaster;
}


// Function to calculate TPR and FPR that returns a Promise-like structure
function calculateTPRandFPRAsync(predictedRaster, processedGTImage, weekFreq, yearFreq, callback) {
  var truePositives = predictedRaster.and(processedGTImage);
  var falsePositives = predictedRaster.and(processedGTImage.not());
  var falseNegatives = processedGTImage.and(predictedRaster.not());
  var trueNegatives = processedGTImage.not().and(predictedRaster.not());

  // Combine all calculations into one reduceRegion call to minimize server requests
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
      var tp = result.TP;
      var fp = result.FP;
      var fn = result.FN;
      var tn = result.TN;
      
      var tpr = tp / (tp + fn);
      var fpr = fp / (fp + tn);
      
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
    // generateROCCurve();
    return;
  }
  
  var combo = parameterCombinations[combinationIndex];
  var weekFreq = combo.weekFreq;
  var yearFreq = combo.yearFreq;
  
  print('Processing combination:', combinationIndex + 1, 'of', totalCombinations, 
        '- WeekFreq:', weekFreq, 'YearFreq:', yearFreq);
  
  var result = processSelectedMask(2023, 5, yearFreq, weekFreq);
  
  // if (result && result.floodResults) {
  //   calculateTPRandFPRAsync(result.floodResults, processedGTImage, weekFreq, yearFreq, function(metrics) {
  //     if (metrics) {
  //       rocResults.push(metrics);
  //       print('Combination', combinationIndex + 1, 'complete - TPR:', metrics.TPR.toFixed(4), 'FPR:', metrics.FPR.toFixed(4));
  //     } else {
  //       print('Failed to calculate metrics for combination', combinationIndex + 1);
  //     }
      
  //     // Process next combination only after current one is complete
  //     processNextCombination(combinationIndex + 1);
  //   });
  // } else {
  //   print('No results for combination', combinationIndex + 1, '- WeekFreq:', weekFreq, 'YearFreq:', yearFreq);
  //   // Continue to next combination immediately if no processing needed
  //   processNextCombination(combinationIndex + 1);
  // }
  
  if(result && result.floodResults){
      processNextCombination(combinationIndex + 1);  
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
