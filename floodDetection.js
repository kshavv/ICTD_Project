//***********************New and improved biweekly flood maps*********************************//


// Defining the area of interest
var districts = ee.FeatureCollection("FAO/GAUL/2015/level2");
var aoi = districts.filter(ee.Filter.and(
  ee.Filter.eq("ADM1_NAME", "Assam"),
  ee.Filter.or(
    ee.Filter.eq("ADM2_NAME", "Kamrup"),
    ee.Filter.eq("ADM2_NAME", "Barpeta"),
    ee.Filter.eq("ADM2_NAME", "Nalbari")
  )
));
Map.centerObject(aoi, 10);


//CONFIGS
var years = ee.List.sequence(2018, 2022);
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

//************************ Find and store the masks for all years *************************************************//
var biWeeklyMasksByYear = years.map(biWeeklyMasks);


//------------------------ Classification starts here -------------------------------

//1. Finding the Perennial Water and non-water pixels
function classifyPerennialAndNonWater(years, monsoonStart, monsoonEnd, aoi, threshold, perennialThreshold) {

  // 2. Calculate the total number of bi-weekly periods across all years.
  var totalBiWeeks = years.map(function(year) {
    var start = ee.Date.fromYMD(year, monsoonStart, 1);
    var end = ee.Date.fromYMD(year, monsoonEnd, 30);
    return end.difference(start, 'week').divide(2).floor();
  }).reduce(ee.Reducer.sum());

  // 3. Sum the bi-weekly water masks to get water presence count.
  var waterPresenceSum = ee.Image(0);
  for (var i = 0; i < years.size().getInfo(); i++) {
    var yearCollection = ee.ImageCollection(biWeeklyMasksByYear.get(i));
    if (yearCollection.size().getInfo() > 0) {
      waterPresenceSum = waterPresenceSum.add(yearCollection.reduce(ee.Reducer.sum()).unmask());
    } else {
      waterPresenceSum = waterPresenceSum.add(ee.Image(0).toByte());
    }
  }

  // 4. Calculate the proportion of time water is present.
  var waterPresenceProportion = waterPresenceSum.divide(ee.Image.constant(totalBiWeeks));
  
  // 5. Classify pixels
  var perennialWater = waterPresenceProportion.gte(perennialThreshold).rename('perennial_water');
  var nonWater = waterPresenceProportion.eq(0).rename('non_water');

  var classification = ee.Image(0)
    .where(perennialWater, 1)
    .where(nonWater, 2);

  return classification.clip(aoi);
}

var perennialAndNonWaterClassification = classifyPerennialAndNonWater(years, monsoonStart, monsoonEnd, aoi, threshold, 0.85);
Map.addLayer(perennialAndNonWaterClassification, { palette: ['000000', '0000FF','00BB00'], min:0, max:2 }, 'Perennial and Non-Water');

// Function to display the selected bi-weekly water mask
// Function to display the selected bi-weekly water mask
// Function to display the selected bi-weekly water mask
function displaySelectedMask() {
  var weekFreq = 0.6; // This is now bi-weekly frequency
  var yearFreq = 0.8;

  var selectedYear = yearSelect.getValue();
  var selectedBiWeek = ee.Number(biWeekSelect.getValue()).toInt();

  // Access the ImageCollection for the selected year
  var selectedYearCollection = ee.ImageCollection(biWeeklyMasksByYear.get(years.indexOf(selectedYear)));

  // Filter the ImageCollection for the selected bi-week
  var selectedMask = selectedYearCollection
    .filterMetadata('biweek', 'equals', selectedBiWeek)
    .first();

  // Check if an image was found BEFORE using it.
  if (selectedMask) {
    // Since we confirmed selectedMask exists, we cast it to an ee.Image
    selectedMask = ee.Image(selectedMask);

    // --- LOGGING CODE RESTORED HERE ---
    var selectedBiWeekStartDate = ee.Date.fromYMD(selectedYear, monsoonStart, 1).advance(selectedBiWeek.multiply(2), 'week');
    print('Processing data for bi-week starting:', selectedBiWeekStartDate.format('YYYY-MM-dd'));
    // ------------------------------------

    var weeklyNonWaterMask = selectedMask.eq(0).rename('weekly_non_water');

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

    var unclassifiedPixels = perennialAndNonWaterClassification.eq(0);
    var newNonWaterPixels = unclassifiedPixels.and(weeklyNonWaterMask);

    var alreadyClassified = perennialAndNonWaterClassification;
    alreadyClassified = alreadyClassified.where(newNonWaterPixels, 2);

    var seasonalWater = waterFrequencyThisBiWeek.gte(weekFreq).and(waterFrequencyThisYear.lt(yearFreq)).rename('seasonal_water').and(alreadyClassified.eq(0));
    var floodWater = waterFrequencyThisBiWeek.lt(weekFreq).and(waterFrequencyThisYear.lt(yearFreq)).rename('flood_water').and(alreadyClassified.eq(0).and(seasonalWater.not()));
    var newPerrenialWater = waterFrequencyThisYear.gte(yearFreq).rename('new_seasonal_water').and(alreadyClassified.eq(0).and(seasonalWater.not().and(floodWater.not())));

    var finalClassification = alreadyClassified.where(seasonalWater, 3).where(floodWater, 4).where(newPerrenialWater, 5);

    Map.layers().set(1, ui.Map.Layer(finalClassification.clip(aoi), {
      palette: ['000000', '0000FF', '00BB00', 'yellow', 'red', '72A1ED'],
      min: 0,
      max: 5
    }, 'Water Classification'));

    Export.image.toDrive({
      image: finalClassification,
      description: 'finalClassification_biweekly',
      fileNamePrefix: 'finalClassification_biweekly',
      region: aoi.geometry(),
      scale: 30,
      crs: 'EPSG:4326',
      maxPixels: 1e13
    });

  } else {
    // This block runs if no image was found for the selection.
    print('No data found for the selected year and bi-week.');
    Map.layers().set(1, ui.Map.Layer(ee.Image().rename('Water Classification'), {}, 'Water Classification', false));
  }
}

// -------------------- UI ----------------------------
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
