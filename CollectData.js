
///////////////////////////////////////////////////////////////
/////////////////PARSING THE DATA//////////////////////////////
///////////////////////////////////////////////////////////////

var parsed = df.map(function(f) {
  var name = ee.String(f.get('Name'));   // e.g. "19W19102022"

  // Extract parts using regex
  var id    = ee.String(name.match('^[0-9]+').get(0));      // "19"
  var type  = ee.String(name.match('W|NW').get(0));         // "W" or "NW"

  // Remove id + type → remaining "19102022"
  var rest = name.replace(id.cat(type), '');

  // Slice DDMMYYYY
  var day   = ee.String(rest.slice(0, 2));   // "19"
  var month = ee.String(rest.slice(2, 4));   // "10"
  var year  = ee.String(rest.slice(4, 8));   // "2022"

  return f.set({
    id: id,
    waterType: type,
    day: day,
    month: month,
    year: year
  });
});

//////////////////////////////////////////////////////////////
//////////////FUNCTION FOR ENRICHING PIXELS //////////////////
//////////////////////////////////////////////////////////////
function getClosestS1Image(targetDate, geom) {

  var start = targetDate.advance(-7, 'day');
  var end   = targetDate.advance(7, 'day');

  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(geom) 
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .map(function(img) {
      var diff = img.date().difference(targetDate, 'day').abs();
      return img.set('timeDiff', diff);
    })
    .sort('timeDiff');

  return ee.Image(s1.first());
}

function getClosestS2Image(targetDate, geom) {

  // Smaller window because S2 has high revisit frequency
  var start = targetDate.advance(-3, 'day');
  var end   = targetDate.advance(3, 'day');

  var s2 = ee.ImageCollection('COPERNICUS/S2_SR')
    .filterBounds(geom) 
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40))
    .map(function(img) {
      var diff = img.date().difference(targetDate, 'day').abs();
      return img.set('timeDiff', diff);
    })
    .sort('timeDiff');

  return ee.Image(s2.first());
}


//////////////////////////////////////////////////////////////
//////////////EXTRACTING PIXEL LEVEL DATA FROM POLYGONS///////
//////////////////////////////////////////////////////////////

//using a lat long based method..(could also use the projection of Sentinel-1)
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

  // ------------------------------------
  // 1. Build exact target date
  // ------------------------------------
  var year  = ee.Number.parse(f.get('year'));
  var month = ee.Number.parse(f.get('month'));
  var day   = ee.Number.parse(f.get('day'));

  var targetDate = ee.Date.fromYMD(year, month, day);
  var geom = f.geometry();

  // ------------------------------------
  // 2. Closest Sentinel-1 (±7 days)
  // ------------------------------------
  var s1 = getClosestS1Image(targetDate, geom);

  var s1DateTime = ee.Date(s1.get('system:time_start'));
  var s1DateStr  = s1DateTime.format('YYYY-MM-dd HH:mm:ss');

  // ------------------------------------
  // 3. Closest Sentinel-2 (±3 days)
  // ------------------------------------
  var s2 = getClosestS2Image(targetDate, geom);

  var s2DateTime = ee.Date(s2.get('system:time_start'));
  var s2DateStr  = s2DateTime.format('YYYY-MM-dd HH:mm:ss');

  // ------------------------------------
  // 4. Stack bands
  // ------------------------------------
  var stacked = s1.select(['VV', 'VH'])
                  .addBands(
                    s2.select(['B2', 'B3', 'B4', 'B8'])
                  );

  // ------------------------------------
  // 5. Sample at pixel location
  // ------------------------------------
  var sampled = stacked.sample({
    region: geom,
    scale: 10
  }).first();

  // ------------------------------------
  // 6. Safe fallback + metadata
  // ------------------------------------
  return ee.Feature(
    ee.Algorithms.If(
      sampled,
      f.copyProperties(sampled, sampled.propertyNames())
       .set({
         s1_datetime_utc: s1DateStr,
         s2_datetime_utc: s2DateStr,
         s1_day_diff: s1.get('timeDiff'),
         s2_day_diff: s2.get('timeDiff')
       }),
      f.set({
        s1_missing: 1,
        s2_missing: 1
      })
    )
  );
});
print(enrichedWithSentinel.first())
// print(enrichedWithSentinel.limit(10))


////////////////////////////////DEBUG////////////////////////////////////////
// Water pixels
var waterPixels = pixelData.filter(
  ee.Filter.eq('waterType', 'W')
);

// Non-water pixels
var nonWaterPixels = pixelData.filter(
  ee.Filter.eq('waterType', 'NW')
);

// Print counts
print('Total water pixels:', waterPixels.size());
print('Total non-water pixels:', nonWaterPixels.size());
print('Total pixels:', pixelData.size());

// Count pixels per polygon id
var pixelCountDict = pixelData.aggregate_histogram('Name');
print('Pixel count per polygon (dict):', pixelCountDict);


var pixelCountFC = ee.FeatureCollection(
  ee.Dictionary(pixelCountDict).keys().map(function(k) {
    return ee.Feature(null, {
      id: k,
      pixel_count: ee.Dictionary(pixelCountDict).get(k)
    });
  })
);

//////////////////////////////////////////////////////////////////////////////


///////////////////////
///sampling new data///
/////////////////////// 

var bufferDist = 20;         // outward buffer (m)
var samplesPerPoly = 20;      // N samples per polygon
var scale = 10;   

// Base image used ONLY for sampling pixel locations
var baseImg = ee.Image.pixelLonLat()
  .addBands(ee.Image.constant(1).rename('dummy'));
  
// Safe empty image for null-guarding Sentinel data
var EMPTY_IMAGE = ee.Image.constant(0).selfMask();

var processPolygon = function(feat) {

  var poly = ee.Feature(feat);
  var geom = poly.geometry();

  // -----------------------------
  // Buffers
  // -----------------------------
  var innerBuffer = geom.buffer(bufferDist);
  var outerBuffer = geom.buffer(bufferDist + 50);
  var ring = outerBuffer.difference(innerBuffer, 1);

  var finalRegion = ee.Geometry(
    ee.Algorithms.If(
      ring.area(1).gt(0),
      ring,
      outerBuffer
    )
  );

  // -----------------------------
  // Sample candidate pixels
  // -----------------------------
  var candidates = baseImg.sample({
    region: finalRegion,
    scale: scale,
    geometries: true
  });

  // -----------------------------
  // Pick N samples
  // -----------------------------
  var selected = candidates
    .randomColumn('rand')
    .sort('rand')
    .limit(samplesPerPoly)
    .map(function(p) {
      return p.set({
        id: poly.get('id'),
        Name: poly.get('Name'),
        waterType: 'NW',
        day: poly.get('day'),
        month: poly.get('month'),
        year: poly.get('year')
      });
    });

  // -----------------------------
  // Sentinel enrichment (NULL-SAFE)
  // -----------------------------
  var enriched = selected.map(function(fp) {

    var targetDate = ee.Date.fromYMD(
      ee.Number.parse(fp.get('year')),
      ee.Number.parse(fp.get('month')),
      ee.Number.parse(fp.get('day'))
    );

    var s1Obj = getClosestS1Image(targetDate, fp.geometry());
    var s2Obj = getClosestS2Image(targetDate, fp.geometry());

    var s1Image = ee.Image(
      ee.Algorithms.If(s1Obj.hasImage, s1Obj.image, EMPTY_IMAGE)
    );

    var s2Image = ee.Image(
      ee.Algorithms.If(s2Obj.hasImage, s2Obj.image, EMPTY_IMAGE)
    );

    var stacked = s1Image.select(['VV', 'VH'])
      .addBands(
        s2Image.select(['B2', 'B3', 'B4', 'B8'])
      );

    var sampled = stacked.sample({
      region: fp.geometry(),
      scale: 10
    }).first();

    // If no data → sampled = null → properties remain null
    return ee.Feature(
      ee.Algorithms.If(
        sampled,
        fp.copyProperties(sampled, sampled.propertyNames()),
        fp
      )
    );
  });

  return enriched;
};



//select only water polys
var waterPolys = parsed.filter(
  ee.Filter.eq('waterType', 'W')
);

var newNWPoints = ee.FeatureCollection(
  waterPolys.map(processPolygon).flatten()
);

print('New NW samples:', newNWPoints.size());
print(newNWPoints.limit(10));


//////////////////////////
//combining the result////
//////////////////////////

var combined = enrichedWithSentinel.merge(newNWPoints);


// ----------------------------------
// Add latitude / longitude
// ----------------------------------
var addLatLon = function(f) {
  var coords = f.geometry().coordinates();
  return f.set({
    longitude: coords.get(0),
    latitude: coords.get(1)
  });
};

var withLatLon = combined.map(addLatLon);

// ----------------------------------
// Fields to keep
// ----------------------------------
var keepList = [
  'id',
  'B2', 'B3', 'B4', 'B8',
  'VH', 'VV',
  'latitude', 'longitude',
  'day', 'month', 'year',
  'waterType'
];

// ----------------------------------
// Safe select (prevents missing-column issues)
// ----------------------------------
var cleaned = withLatLon.map(function(f) {
  return ee.Feature(null).set(
    ee.Dictionary.fromLists(
      keepList,
      keepList.map(function(k) {
        return f.get(k);
      })
    )
  );
});


// ----------------------------------
// Remove geometry (CSV-friendly)
// ----------------------------------
var cleanedNoGeom = cleaned.map(function(f) {
  return f.setGeometry(null);
});


// ----------------------------------
// Final checks
// ----------------------------------
print('Final export rows:', cleanedNoGeom.size());
print('Example row:', cleanedNoGeom.first());


// ----------------------------------
// EXPORT
// ----------------------------------
// Export.table.toDrive({
//   collection: cleanedNoGeom,
//   description: 'final_pixel_dataset_cleaned',
//   folder: 'keshav_sparsh/Tiffs',
//   fileNamePrefix: 'pixel_data_v1',
//   fileFormat: 'CSV'
// });























///////////////////////////////////////////////////////
/////////////DEBUG POLYGON LEVEL///////////////////////
///////////////////////////////////////////////////////
var list = parsed.toList(parsed.size());
// Example: print feature with index 3 (the 4th feature)


// ---- Count Water (W) and Non-water (N) ----
var water     = parsed.filter(ee.Filter.eq('waterType', 'W'));
var nonWater  = parsed.filter(ee.Filter.eq('waterType', 'NW'));
var others    = parsed.filter(
                  ee.Filter.neq('waterType', 'W')
                ).filter(
                  ee.Filter.neq('waterType', 'NW')
                );
// Print results
print('Water count:', water.size());
print('Non-water count:', nonWater.size());
print('Other types count:', others.size());   
print('Total polygons:', parsed.size());

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
