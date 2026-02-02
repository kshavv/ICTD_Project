
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
  
  var area_m2 = f.geometry().area(1);
  
  return f.set({
    id: id,
    waterType: type,
    day: day,
    month: month,
    year: year,
    poly_area_m2: area_m2
  });
});


////////////////////////////////////////////////////////////////////////////////
//////////////FUNCTION FOR ENRICHING PIXELS WITH SENTINEL DATA//////////////////
////////////////////////////////////////////////////////////////////////////////
function getClosestS1Image(targetDate, geom) {

  // ±7 day search window
  var start = targetDate.advance(-7, 'day');
  var end   = targetDate.advance(7, 'day');

  var s1col = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(geom)
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains(
      'transmitterReceiverPolarisation', 'VV'
    ))
    .filter(ee.Filter.listContains(
      'transmitterReceiverPolarisation', 'VH'
    ))
    .map(function(img) {
      var diff = img.date()
        .difference(targetDate, 'day')
        .abs();
      return img
        .select(['VV','VH'])          // select HERE
        .set('timeDiff', diff);
    })
    .sort('timeDiff');

  // Fallback WITH SAME BAND NAMES
  var empty = ee.Image.constant([0,0])
    .rename(['VV','VH'])
    .selfMask();

  return ee.Image(
    ee.Algorithms.If(
      s1col.size().gt(0),
      s1col.first(),
      empty
    )
  );
}

function getClosestS2Image(targetDate, geom) {

  var start = targetDate.advance(-3, 'day');
  var end   = targetDate.advance(3, 'day');

  var s2col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(geom)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
    .map(function(img) {
      var diff = img.date()
        .difference(targetDate, 'day')
        .abs();
      return img
        .select(['B2','B3','B4','B8','B5','B6','B7','B8A'])   // select HERE
        .set('timeDiff', diff);
    })
    .sort('timeDiff');

  // Fallback WITH SAME BAND NAMES
  var empty = ee.Image.constant([0,0,0,0,0,0,0,0])
    .rename(['B2','B3','B4','B8','B5','B6','B7','B8A'])
    .selfMask();

  return ee.Image(
    ee.Algorithms.If(
      s2col.size().gt(0),
      s2col.first(),
      empty
    )
  );
}

///////////////////////////////////////////////////////////////
//////////////CONVERTING POLYGONS INTO INDIVIDUAL PIXELS///////
///////////////////////////////////////////////////////////////
//Using a Lat Long based method..(could also use the projection of Sentinel-1)
var base = ee.Image.pixelLonLat().rename(['longitude', 'latitude']).addBands(ee.Image.constant(1).rename('maskBase'));
var fc = parsed.map(function(f) {
  return f.set('poly_id', f.get('Name'));
});



var raster = ee.Image().byte().paint({
  featureCollection: fc,
  color: 1
}).rename('mask');


var perPolyLimited = ee.FeatureCollection(
  parsed.map(function(f) {

    var geom = f.geometry();

    var pixels = base.sample({
      region: geom,
      scale: 10,
      geometries: true
    });

    return pixels
      .randomColumn('rand')
      .sort('rand')
      .limit(10)
      .map(function(p) {
        return p.set({
          id: f.get('id'),
          Name: f.get('Name'),
          day: f.get('day'),
          month: f.get('month'),
          year: f.get('year'),
          waterType: f.get('waterType'),
          poly_area_m2: f.get('poly_area_m2') 
        });
      });
  })
).flatten();

// print('Histogram by polygon Name:',
//   perPolyLimited.aggregate_histogram('poly_id')
// );

print(perPolyLimited.limit(20));


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////








////////////////////////////////////////
///SAMPLING NEW NON WATER DATA POINTS///
/////////////////////// ////////////////

var bufferDist = 30;         // outward buffer (m)
var samplesPerPoly = 10;      // N samples per polygon
var scale = 10;   

// Base image used ONLY for sampling pixel locations
var baseImg = ee.Image.pixelLonLat()
  .addBands(ee.Image.constant(1).rename('dummy'));
  
var processPolygon = function(feat) {

  var poly = ee.Feature(feat);
  var geom = poly.geometry();

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

  var candidates = baseImg.sample({
    region: finalRegion,
    scale: scale,
    geometries: true
  });

  return candidates
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
        year: poly.get('year'),
        poly_area_m2: poly.get('poly_area_m2')
      });
    });
};


//select only the water polys
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

var combined = perPolyLimited.merge(newNWPoints);

// print('combined data points')
// print(combined.limit(10));


////////////////////////////////
////add sentinel information//// 
////////////////////////////////
var finalDataset = combined.map(function(f) {

  // ------------------------------------
  // 1. Build exact target date
  // ------------------------------------
  var year  = ee.Number.parse(f.get('year'));
  var month = ee.Number.parse(f.get('month'));
  var day   = ee.Number.parse(f.get('day'));

  var targetDate = ee.Date.fromYMD(year, month, day);
  var geom = f.geometry();

  // -------------------------------------------------------------------
  // 2. Closest Sentinel-1 (±7 days) and Closest Sentinel-2 (±3 days)
  // -------------------------------------------------------------------
  var s1 = getClosestS1Image(targetDate, geom);

  var s1DateTime = ee.Date(s1.get('system:time_start'));
  var s1DateStr  = s1DateTime.format('YYYY-MM-dd HH:mm:ss');
  
  var s2 = getClosestS2Image(targetDate, geom);
  var s2DateTime = ee.Date(s2.get('system:time_start'));
  var s2DateStr  = s2DateTime.format('YYYY-MM-dd HH:mm:ss');
  
  // check if the images exist for the given period
  var s1Exists = s1.bandNames().size().gt(0);
  var s2Exists = s2.bandNames().size().gt(0);

  // ------------------------------------
  // 3. Stacking  the bands if both image exists for the given period
  // ------------------------------------
    return ee.Feature(
    ee.Algorithms.If(
      s1Exists.and(s2Exists),

      (function () {
        var s1DateTime = ee.Date(s1.get('system:time_start'));
        var s2DateTime = ee.Date(s2.get('system:time_start'));

        var stacked = s1.select(['VV', 'VH'])
                         .addBands(
                           s2.select(['B2','B3','B4','B8','B5','B6','B7','B8A'])
                         );

        var sampled = stacked.sample({
          region: geom,
          scale: 10
        }).first();

        return ee.Algorithms.If(
          sampled,
          f.copyProperties(sampled, sampled.propertyNames())
           .set({
             s1_datetime_utc: s1DateTime.format('YYYY-MM-dd HH:mm:ss'),
             s2_datetime_utc: s2DateTime.format('YYYY-MM-dd HH:mm:ss'),
             s1_day_diff: s1.get('timeDiff'),
             s2_day_diff: s2.get('timeDiff'),
             sample_missing:0
           }),
          f.set({ sample_missing: 1 })
        );
      })(),
    
      //if condition not passed(either of data not present)
      f.set({
        s1_missing: s1Exists.not(),
        s2_missing: s2Exists.not()
      })
    )
  );
});

// print(finalDataset.first())
print('Enriched With Sentinel points :', finalDataset.size());
print(finalDataset.limit(20))













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

var withLatLon = finalDataset.map(addLatLon);

// ----------------------------------
// Fields to keep
// ----------------------------------
var keepList = [
  'id','Name',
  'B2', 'B3', 'B4', 'B8','B5','B6','B7','B8A',
  'VH', 'VV',
  'latitude', 'longitude',
  'day', 'month', 'year',
  'poly_area_m2',        
  'waterType','s1_day_diff','s2_day_diff','sample_missing',
  's1_datetime_utc','s2_datetime_utc'
];


print(
  'Area consistency check:',
  combined
    .aggregate_histogram('poly_area_m2')
);
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
Export.table.toDrive({
  collection: cleanedNoGeom,
  description: 'final_pixel_dataset_cleaned',
  folder: 'keshav_sparsh/Tiffs',
  fileNamePrefix: 'pixel_data_v2',
  fileFormat: 'CSV'
});

// Export.table.toDrive({
//   collection: cleanedNoGeom.limit(10),
//   description: 'final_pixel_dataset_cleaned',
//   folder: 'keshav_sparsh\\Tiffs',
//   fileNamePrefix: 'pixel_data_v1_test',
//   fileFormat: 'CSV'
// });

























////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



////////////////////////////////DEBUG////////////////////////////////////////
print('-------------Start of Debug Block-------------')
var polyName = '101W05042025';
// 08W05102022
//04W29102019
var poly = ee.Feature(
  parsed.filter(ee.Filter.eq('Name', polyName)).first()
);

print('Selected polygon:', poly);

Map.centerObject(poly, 16);

var geom = poly.geometry();

var innerBuffer = geom.buffer(bufferDist);
var outerBuffer = geom.buffer(bufferDist + 50);

var ring = ee.Geometry(
  ee.Algorithms.If(
    outerBuffer.difference(innerBuffer, 1).area(1).gt(0),
    outerBuffer.difference(innerBuffer, 1),
    outerBuffer
  )
);

var baseImg = ee.Image.pixelLonLat()
  .addBands(ee.Image.constant(1).rename('dummy'));

var candidates = baseImg.sample({
  region: ring,
  scale: scale,
  geometries: true
});

var nwPoints = candidates
  .randomColumn('rand')
  .sort('rand')
  .limit(samplesPerPoly);


//polygon
Map.addLayer(poly, {
  color: 'blue'
}, 'Original Polygon');

//Inner buffer
Map.addLayer(innerBuffer, {
  color: 'green'
}, 'Inner Buffer'); 

//outer buffer
Map.addLayer(outerBuffer, {
  color: 'red'
}, 'Outer Buffer');

// ring area
Map.addLayer(ring, {
  color: 'orange'
}, 'Sampling Ring');


//NW samples
Map.addLayer(nwPoints, {
  color: 'black',
  pointSize: 4
}, 'NW Sampled Points');



var fc = ee.FeatureCollection([poly]);

var base = ee.Image.pixelLonLat()
  .rename(['lon', 'lat'])
  .addBands(ee.Image.constant(1).rename('maskBase'));
  
  var rasterMask = ee.Image()
  .byte()
  .paint({
    featureCollection: fc,
    color: 1
  })
  .rename('mask');
  
  
  var pixelData = base
  .addBands(rasterMask)
  .sampleRegions({
    collection: fc,
    scale: 10,
    geometries: true
  });


Map.addLayer(poly.geometry(), { color: 'blue' }, 'Polygon');

Map.addLayer(pixelData, {
  color: 'cyan',
  pointSize: 2
}, 'All Pixels Inside Polygon');

print('Pixel count inside polygon:', pixelData.size());
print(pixelData.limit(10));
//////////////////////////////////////////////////////////////////////////////
