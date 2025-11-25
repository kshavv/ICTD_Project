
var parsed = df.map(function(f) {
  var name = ee.String(f.get('Name'));   // e.g., "1W03/22"
  var id = ee.String(name.match('^[0-9]+').get(0));  
  var wtype = ee.String(name.match('W|NW').get(0));
  var after = name.replace(id.cat(wtype), '');  // e.g. "04/22"
  
  var parts = after.split('/');
  var month = ee.String(parts.get(0));   // "04"
  var year  = ee.String(parts.get(1));   // "22"

  return f.set({
    id: id,
    waterType: wtype,
    month: month,
    year: year
  });
});


function getS1Image(start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VV', 'VH'])
    .median();
}

function getS2Image(start, end) {
  return ee.ImageCollection('COPERNICUS/S2_SR')
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40))
    .select(['B2', 'B3', 'B4', 'B8'])
    .median();
}



var list = parsed.toList(parsed.size());
// Example: print feature with index 3 (the 4th feature)
print(list.get(1));

var first = parsed.first();
var firstFC = ee.FeatureCollection([first]);

// ---- Count Water (W) and Non-water (N) ----
var water     = parsed.filter(ee.Filter.eq('waterType', 'W'));
var nonWater  = parsed.filter(ee.Filter.eq('waterType', 'NW'));
var others    = parsed.filter(
                  ee.Filter.neq('waterType', 'W')
                ).filter(
                  ee.Filter.neq('waterType', 'N')
                );
// Print results
print('Water count:', water.size());
print('Non-water count:', nonWater.size());
print('Other types count:', others.size());   
print('Total polygons:', parsed.size());


var base = ee.Image.pixelLonLat().rename(['lon', 'lat'])
             .addBands(ee.Image.constant(1).rename('maskBase'));

var fc = parsed;    // your polygon FC

var raster = ee.Image().byte().paint({
  featureCollection: fc,
  color: 1
}).rename('mask');

var pixelData = base.addBands(raster)
  .sampleRegions({
    collection: fc,
    scale: 10,
    geometries: true
});


////////////////////////////////
////add sentinel information//// 
////////////////////////////////
var enrichedWithSentinel = pixelData.map(function(f) {

  var year  = ee.Number.parse(f.get('year')).add(2000); 
  var month = ee.Number.parse(f.get('month'));          

  var start = ee.Date.fromYMD(year, month, 1);
  var end   = start.advance(1, 'month');

  // Monthly median S1/S2
  var s1 = getS1Image(start, end);
  var s2 = getS2Image(start, end);

  var stacked = s1.addBands(s2);

  // Sample the S1/S2 values at the pixel geometry
  var sampled = stacked.sample({
    region: f.geometry(),
    scale: 10
  }).first();

  // If no S1/S2 image exists → sampled = null
  // Handle safely:
  var result = ee.Algorithms.If(
    sampled,
    f.copyProperties(sampled, sampled.propertyNames()),
    f // fallback: return original pixel without S1/S2
  );

  return ee.Feature(result);
});

// Inspect output
// Total pixels
print('Total pixel count:', enrichedWithSentinel.size());

// Water pixels (from water polygons)
var waterPixels = enrichedWithSentinel.filter(ee.Filter.eq('waterType', 'W'));
print('Water pixel count:', waterPixels.size());

// Non-water pixels
var nonWaterPixels = enrichedWithSentinel.filter(ee.Filter.eq('waterType', 'NW'));
print('Non-water pixel count:', nonWaterPixels.size());
print(enrichedWithSentinel.limit(10));


// Export.table.toDrive({
//   collection: pixelData,
//   description: 'pixel_table_export',
//   fileFormat: 'CSV'
// });


//vis
// Create a palette (up to ~20 colors)
Map.addLayer(raster.selfMask(), {palette: ['yellow']}, 'All polygon pixels');
Map.centerObject(fc, 12);


var styled = firstFC.style({
  color: 'FF0000',
  fillColor: 'FF000055',
  width: 3
});

Map.addLayer(styled, {}, 'First Polygon (styled)');
Map.centerObject(first, 20);
