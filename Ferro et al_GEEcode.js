/*
  The GEE code imports a few assets, see here:
  https://code.earthengine.google.com/c89c63ebd4460e2f527f71cb44d6b807
*/

// 1. Get data information

// 1a. Upload region of interest


// 1b. Get SPEI dataset
var speiChaco = ee.ImageCollection("CSIC/SPEI/2_8")
    .filterDate('2000-01-01', '2023-01-01')
    .map(function(img){
 return img.select('SPEI_06_month') //el muestro original se hizo con SPEI_12
           .rename('spei')
           .clip(hexChaco)
           .set({'system:time_start': ee.Date(img.get('system:time_start'))})});

// Set the visualization ranges and color palette.
var visParams = {
  min: -2.33,
  max:  2.33,
  palette: [
    '8b1a1a', 'de2929', 'f3641d',
    'fdc404', '9afa94', '03f2fd',
    '12adf3', '1771de', '00008b',
  ]
};

// List of years under analysis.
var years = ee.List.sequence(2000, 2022); //sequence of numbers from start to end (inclusive)
var year_id = ee.List.sequence(0, years.size().subtract(1))

// Get humid SPEI season.
var spei_anual_humedo = ee.ImageCollection(years.map(function(y){
  var startDate = ee.Date.fromYMD(ee.Number(y).add(-1), 11, 01),
      endDate = ee.Date.fromYMD(ee.Number(y), 04, 01),
      img = speiChaco.filterDate(startDate, endDate)
                     .mean()
                     .set("year", y);
  return img
}))
print(spei_anual_humedo, 'spei_anual_humedo')
Map.addLayer(spei_anual_humedo, {}, 'spei_anual_humedo', false)


// 1c. Map Biomas Chaco (www.chaco.mapbiomas.org)
var mb_roi = mb.clip(roi)

var bands = mb_roi.bandNames()
var list = bands.map(function(n) { return mb_roi.select([n]) })
var collection = ee.ImageCollection.fromImages(list)

var escalatemp_list = collection.toList(collection.size())

// Get elements from the list (15 to 38) which correspond to the study period.
var escalaTemp = escalatemp_list.slice(15, 38); //start index, inclusive, and end index, exclusive

var col = ee.ImageCollection(escalaTemp)

// Reclass.
var fun_reclas = function(img){
  return img.expression('b(0) <= 3 ? 1' +
                        ': b(0) == 4 ? 1' +
                        ': b(0) == 45 ? 1' +
                        ': b(0) == 42 ? 4' + ': b(0) == 43 ? 4' + ': b(0) == 44 ? 4' + ': b(0) == 11 ? 4' +
                        ': b(0) == 15 ? 5' +
                        ': b(0) == 19 ? 6' +  ': b(0) == 57 ? 6' +  ': b(0) == 58 ? 6' + ': b(0) == 36 ? 6' +
                        ': 0')
            .rename('lulc')
            .clip(roi)
};
// (1) Forest = Closed, open, and scattered woody vegetation
// (4) Grassland
// (5) Pastures
// (6) Agriculture (includes shrub crops)
// (0) Other = All pixels not included in these categories (floodable woody vegetation, forest plantations)

var reclas = col.map(fun_reclas);

Map.addLayer(reclas, {}, 'reclas', false)

// Count pixels for each land cover class.
var resample = reclas.map(function(img){
  return img.reduceRegions({collection: roi,
                            reducer: ee.Reducer.frequencyHistogram(). unweighted(),
                            scale: 30})
}).flatten();


var propList = function(feat){
var dict = ee.Dictionary(feat.get('histogram'));
var llaves = dict.keys()
  return feat.set('llaves', llaves)
}

var result = resample.map(propList)

var coberturas = result.map(function(feat){
  var llaves = ee.List(feat.get('llaves'));
  var dict = ee.Dictionary(feat.get('histogram'));
  return feat.set({Bosque: ee.Algorithms.If(llaves.contains('1'), ee.Number(dict.get('1')).int(), 0),
                   Pastizal:ee.Algorithms.If(llaves.contains('4'), ee.Number(dict.get('4')).int(), 0),
                   Pastura:ee.Algorithms.If(llaves.contains('5'), ee.Number(dict.get('5')).int(), 0),
                   Agricultura:ee.Algorithms.If(llaves.contains('6'), ee.Number(dict.get('6')).int(), 0),
                   Otro:ee.Algorithms.If(llaves.contains('0'), ee.Number(dict.get('0')).int(), 0)
  });
});

// Assign year.
var cob_lista = reclas.toList(reclas.size())

var cob_anual = ee.ImageCollection(year_id.map(function(y){
  var img = ee.Image(ee.List(cob_lista).get(y)).set("year", years.get(y));
  return img
}))
print(cob_anual, 'cob_anual')


// 1d. Burned Area (MODIS-MCD64A1.061)
var bands = colFuegos.bandNames()
var list = bands.map(function(n) { return colFuegos.select([n]).rename('fuegos') })
var fuegos = ee.ImageCollection.fromImages(list)

// Assign year.
var fuegos_lista = fuegos.toList(fuegos.size())

var fuegos_anual = ee.ImageCollection(year_id.map(function(y){
  var img = ee.Image(ee.List(fuegos_lista).get(y)).set("year", years.get(y));
  return img
}))



// 2. Resample
// Minimun distance and number of points per year.
var min_distance = ee.Number(5000)
var n_pix_year = ee.Number(2000)

var year2 = ee.List.sequence(2001, 2022); //End is inclusive.

var datos = ee.FeatureCollection(year2.map(function(y){
  y = ee.Number(y);

  // Multiband image to sample points.
  var img = fuegos_anual.filterMetadata("year", "equals", y).first()
                  .addBands([
                     ee.Image.constant(ee.Number(y)),
                      cob_anual.filterMetadata("year", "equals", y).first(),
                      spei_anual_humedo.filterMetadata("year", "equals", y).first(),
                      spei_anual_seco.filterMetadata("year", "equals", y).first(),
                   //   fwi_anual.filterMetadata("year", "equals", y).first(),
                      shannon_anual.filterMetadata("year", "equals", y).first(),
                      cob_anual.filterMetadata("year", "equals", y.add(-1)). first(),
                  //    cob_anual.filterMetadata("year", "equals", y.add(1)). first(),
                      fuegos_anual.filterMetadata("year", "equals", y.add(-1)). first()
                  //    fuegos_anual.filterMetadata("year", "equals", y.add(1)). first()
                    ]).rename(["fuego_focal", "anio", "cob_focal", "speiHumedo", "speiSeco", "shannon", "cob_prev", "fuego_prev"])

    // Generation of random points with a minimun distance.
    // Make points.
    var seed = ee.Number(y),
    proj = ee.Projection("EPSG:5346").atScale(min_distance), //proyeccion plana para medir distancia entre puntos

    // Use proyection posgar 2007.
    cells = ee.Image.random(seed).multiply(100000).int()
              .clip(roi).reproject(proj);

// The random seed is a key factor in the generation of random points. If the seed remains
//constant between runs, then the sequence of pseudorandom numbers generated will be identical,
//resulting in the same point patterns. Furthermore, if all other parameters and processes are also kept
//constant, the result of point generation will be consistent across runs.

    // Generate another random image and select the maximum random value
    // in each grid cell as the sample point.
    var random = ee.Image.random(seed).multiply(1000000).int();

    // Create a mask to remove every pixel with even coordinates in either X or Y.
    // Using the double not to avoid floating point comparison issues.
    var mask = ee.Image.pixelCoordinates(proj)
                 .expression("!((b('x') + 0.5) % 2 != 0 || (b('y') + 0.5) % 2 != 0)"),
        strictCells = cells.updateMask(mask).reproject(proj),

    // Get maximum pixels within each larger cell.
    strictMax = strictCells.addBands(random)
                           .reduceConnectedComponents(ee.Reducer.max());

    // Find all the points that are local maximums.
    var strictPoints = random.eq(strictMax).selfMask()
                             .clip(roi);

    // Turn into points
    // The SPEI pixel is resampled to 500m, using interpolation based on the nearest centroids.
    var strictSample = img.updateMask(strictPoints.mask()).sampleRegions({
      collection: roi,
      scale: 500,
      geometries: true
    });

  return strictSample
}));


// 3. Export data.