// ============================================================
//  Snow Depth Estimation — SAR + Optical Semi-Empirical Model
//  Sensor   : Sentinel-1 GRD (SAR) + Sentinel-2 SR (Optical)
//  Method   : VH/VV backscatter ratio × vegetation-corrected
//             empirical model, restricted to NDSI snow mask
//  Author   : Shahid Shuja Shafai  <shahidshafai@gmail.com>
//  Lab      : Himalayan Cryospheric Research Lab,
//             University of Kashmir
// ============================================================

// ── DESCRIPTION ─────────────────────────────────────────────
//
// Snow depth (SD) estimation from spaceborne SAR is based on
// the physical relationship between snowpack volume scattering
// and cross-polarisation backscatter. As snow depth increases,
// volume scattering within the snowpack raises VH relative to VV,
// making the VH/VV ratio a proxy for snow water content and depth.
//
// This script implements a semi-empirical SD model:
//
//   SD (cm) = [100 × SI_cumulative × a] / [1 − b × (NDVI − MSAVI)]
//
// Where:
//   SI_cumulative — cumulative sum of VH/VV ratios across five
//     Sentinel-1 acquisitions spanning the snow season. Cumulative
//     accumulation integrates the seasonal backscatter signal rather
//     than relying on a single date, improving robustness against
//     temporal noise and wet-snow effects on individual scenes.
//
//   a = 0.82 — empirical scaling coefficient
//   b = 0.6  — vegetation dampening factor
//
//   (NDVI − MSAVI) — vegetation correction sub-index. MSAVI
//     suppresses the soil background signal, so the residual
//     (NDVI − MSAVI) isolates the vegetation canopy contribution
//     not accounted for by MSAVI alone. This term corrects the
//     denominator for forest/shrub canopy effects on the radar
//     signal, preventing overestimation of SD under vegetation.
//     Alpine vegetation in the study region shows negligible
//     phenological change year-round, so optical acquisitions
//     from any cloud-free window within the year are valid for
//     deriving a stable vegetation correction layer.
//
// PRE-PROCESSING:
//   • Terrain correction (volumetric model) applied to all S1
//     images to remove topographic distortion in mountainous terrain
//   • Speckle filtering via 100 m focal median before ratio computation
//   • Snow spatial mask derived from optical NDSI > 0.4 restricts
//     SD estimation to snow-covered pixels only — the model is not
//     physically meaningful over bare soil, rock, or vegetation
//
// OUTPUT: Snow depth map (cm) clipped to optical snow mask,
//         exportable as cloud-optimised GeoTIFF.
//
// ── HOW TO USE ───────────────────────────────────────────────
//  1. Replace AOI with your study area (see options below)
//  2. Adjust date windows for your target snow season
//  3. Adjust S1 orbit numbers for your region
//  4. Run — check Console for collection sizes and SD statistics
//  5. Export via Tasks tab
//
// ─────────────────────────────────────────────────────────────


// ============================================================
//  ★  USER-CONFIGURABLE PARAMETERS — EDIT HERE  ★
// ============================================================

// ── AOI ──────────────────────────────────────────────────────
// Replace with your own study area. Options:
//   a) Drawn geometry:  var AOI = geometry;
//   b) Uploaded asset:  var AOI = ee.FeatureCollection('users/you/boundary');
//   c) Admin boundary:  var AOI = ee.FeatureCollection('FAO/GAUL/2015/level2')
//                                   .filter(ee.Filter.eq('ADM2_NAME','YourDistrict'));
var AOI = table2;   // <-- REPLACE with your study area

// ── SENTINEL-2 SNOW SEASON WINDOW ────────────────────────────
// Used for: NDSI snow mask + NDVI/MSAVI vegetation correction
// Alpine vegetation is phenologically stable year-round, so any
// cloud-free winter window is valid for the vegetation correction.
var S2_START = '2021-12-01';
var S2_END   = '2022-03-30';
var S2_MAX_CLOUD = 15;   // % — lower is better for snow mapping

// ── SENTINEL-1 DATE WINDOW ───────────────────────────────────
// Should span the target snow accumulation season.
// Ensure at least 5 acquisitions exist (check Console output).
var S1_START = '2021-01-01';
var S1_END   = '2021-03-15';

// ── S1 ORBIT NUMBER ──────────────────────────────────────────
// Relative orbit covering your AOI. Check S1 orbit finder tool.
var S1_ORBIT = 34;

// ── NDSI SNOW MASK THRESHOLD ─────────────────────────────────
// Pixels with NDSI > threshold are classified as snow and used
// to spatially restrict the SD model. 0.4 is conservative.
var NDSI_THRESHOLD = 0.4;

// ── SD MODEL COEFFICIENTS ────────────────────────────────────
// Empirically determined. Adjust for regional calibration.
var SD_COEFF_A = 0.82;   // scaling coefficient
var SD_COEFF_B = 0.6;    // vegetation dampening factor

// ── EXPORT SETTINGS ──────────────────────────────────────────
var EXPORT_FOLDER      = 'GEE_Snow_Depth_Outputs';
var EXPORT_DESCRIPTION = 'SnowDepth_SAR_Optical';
var EXPORT_FILENAME    = 'SnowDepth_' + S1_START.slice(0, 7);
var EXPORT_SCALE       = 20;   // metres

// ============================================================
//  END OF USER PARAMETERS
// ============================================================


// ── Map Initialisation ───────────────────────────────────────
Map.centerObject(AOI, 10);
Map.setOptions('SATELLITE');

// Shared colour palette (cold→hot, 30 stops)
var SD_PALETTE = [
  '#000080','#0000a7','#0000cf','#0000f7','#000dff',
  '#0030ff','#0054ff','#0077ff','#009aff','#00bdff',
  '#00e0fb','#18ffdf','#34ffc2','#51ffa6','#6dff8a',
  '#8aff6d','#a6ff51','#c2ff34','#dfff18','#fbf100',
  '#ffd000','#ffb000','#ff8f00','#ff6e00','#ff4e00',
  '#ff2d00','#f70d00','#cf0000','#a70000','#800000'
];


// ── Terrain Correction — Sentinel-1 ─────────────────────────
// Volumetric model corrects for topographic distortion in
// mountainous terrain. Applies layover and shadow masking.
// Reference: Vollrath et al. 2020, Remote Sensing.
function terrainCorrection(image) {
  var imgGeom    = image.geometry();
  var srtm       = ee.Image('USGS/SRTMGL1_003').clip(imgGeom);
  var sigma0Pow  = ee.Image.constant(10).pow(image.divide(10.0));

  var theta_i = image.select('angle');
  var phi_i   = ee.Terrain.aspect(theta_i)
    .reduceRegion(ee.Reducer.mean(), theta_i.get('system:footprint'), 1000)
    .get('aspect');

  var alpha_s = ee.Terrain.slope(srtm).select('slope');
  var phi_s   = ee.Terrain.aspect(srtm).select('aspect');
  var phi_r   = ee.Image.constant(phi_i).subtract(phi_s);

  var phi_rRad   = phi_r.multiply(Math.PI / 180);
  var alpha_sRad = alpha_s.multiply(Math.PI / 180);
  var theta_iRad = theta_i.multiply(Math.PI / 180);
  var ninetyRad  = ee.Image.constant(90).multiply(Math.PI / 180);

  // Slope steepness in range (eq. 2) and azimuth (eq. 3)
  var alpha_r  = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan();
  var alpha_az = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan();

  // Local incidence angle (eq. 4)
  var theta_lia    = (alpha_az.cos().multiply((theta_iRad.subtract(alpha_r)).cos())).acos();
  var theta_liaDeg = theta_lia.multiply(180 / Math.PI);

  // Gamma-naught flat
  var gamma0   = sigma0Pow.divide(theta_iRad.cos());
  var gamma0dB = ee.Image.constant(10).multiply(gamma0.log10());

  // Volumetric model
  var volModel       = ((ninetyRad.subtract(theta_iRad).add(alpha_r)).tan())
                         .divide((ninetyRad.subtract(theta_iRad)).tan()).abs();
  var gamma0_Volume  = gamma0.divide(volModel);
  var gamma0_VolumedB = ee.Image.constant(10).multiply(gamma0_Volume.log10());

  // Layover mask (slope > radar viewing angle) and shadow mask (LIA > 85°)
  var layover = alpha_r.multiply(180 / Math.PI).lt(theta_i);
  var shadow  = theta_liaDeg.lt(85);

  var ratio = gamma0_VolumedB.select('VV').subtract(gamma0_VolumedB.select('VH'));

  var output = gamma0_VolumedB
    .addBands(ratio).addBands(alpha_r).addBands(phi_s)
    .addBands(theta_iRad).addBands(layover).addBands(shadow)
    .addBands(gamma0dB);

  return image.addBands(
    output.select(['VV', 'VH', 'slope_1', 'slope_2'], ['VV', 'VH', 'layover', 'shadow']),
    null, true
  );
}


// ── Speckle Filter ───────────────────────────────────────────
// 100 m focal median applied to all bands before ratio computation
function filterSpeckles(img) {
  return img.addBands(img.focal_median(100, 'square', 'meters'));
}


// ── Cloud Masking — Sentinel-2 SR ────────────────────────────
function maskS2clouds(image) {
  var qa            = image.select('QA60');
  var cloudBitMask  = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var clearMask = qa.bitwiseAnd(cloudBitMask).eq(0)
                    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(clearMask).divide(10000);
}


// ── Sentinel-2 — Snow Season Composite ──────────────────────
print('── Building Sentinel-2 winter composite ──');
print('Window: ' + S2_START + ' → ' + S2_END);

var s2Collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate(S2_START, S2_END)
  .filterBounds(AOI)
  .map(function(image) { return image.clip(AOI); })
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', S2_MAX_CLOUD))
  .map(maskS2clouds)
  .select(['B2', 'B3', 'B4', 'B8', 'B11']);

print('S2 collection size:', s2Collection.size());

var s2 = s2Collection.median().clip(AOI);

Map.addLayer(s2, { min: 0.05, max: 0.35, bands: ['B4', 'B3', 'B2'] }, 'S2 True Colour');
Map.addLayer(s2, { min: 0.05, max: 0.76, bands: ['B8', 'B4', 'B3'] }, 'S2 NIR False Colour');


// ── NDSI Snow Mask ───────────────────────────────────────────
// NDSI = (Green − SWIR1) / (Green + SWIR1)
// Pixels > NDSI_THRESHOLD classified as snow and used to
// spatially restrict SD model to physically meaningful areas.
var ndsi = s2.normalizedDifference(['B3', 'B11']).rename('NDSI');
Map.addLayer(ndsi, { min: -0.3, max: 0.8 }, 'NDSI');

var snowBinary = ndsi.gt(NDSI_THRESHOLD).clip(AOI);
var snowMask   = ndsi.updateMask(snowBinary).toInt().clip(AOI);

// Vectorise snow mask to define spatial extent for SD clipping and export
var snowPoly = snowMask.reduceToVectors({
  reducer  : ee.Reducer.countEvery(),
  geometry : AOI,
  scale    : 10,
  maxPixels: 1e13
}).geometry();

Map.addLayer(snowMask, {}, 'NDSI Snow Mask (threshold > ' + NDSI_THRESHOLD + ')');


// ── Vegetation Correction Indices ────────────────────────────
// NDVI and MSAVI computed from the same winter composite.
// Alpine vegetation is phenologically stable year-round, making
// cloud-free winter imagery valid for vegetation correction.
//
// MSAVI: Modified Soil-Adjusted Vegetation Index
//   Reduces soil background effects compared to NDVI.
//   Note: GEE expressions do not support pow() — using explicit x*x form.
var ndvi = s2.normalizedDifference(['B8', 'B4']).rename('NDVI').clamp(0.2, 0.8);

var msavi = s2.expression(
  '(2 * NIR + 1 - sqrt((2 * NIR + 1) * (2 * NIR + 1) - 8 * (NIR - RED))) / 2',
  { NIR: s2.select('B8'), RED: s2.select('B4') }
).rename('MSAVI').clamp(0.2, 0.8);

// Vegetation sub-index: residual canopy signal not captured by MSAVI
// Used to correct SD denominator for vegetation effects
var subIndex = ndvi.subtract(msavi);


// ── Sentinel-1 — VH/VV Ratio Stack ──────────────────────────
print('── Building Sentinel-1 SAR stack ──');
print('Window: ' + S1_START + ' → ' + S1_END);

var s1Collection = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
  .filterMetadata('resolution_meters', 'equals', 10)
  .filterBounds(AOI)
  .filterDate(S1_START, S1_END)
  .filter(ee.Filter.inList('relativeOrbitNumber_start', [S1_ORBIT]))
  .map(terrainCorrection)
  .map(filterSpeckles)
  .select(['VV', 'VH', 'angle', 'layover', 'shadow']);

print('S1 collection size (must be ≥ 5):', s1Collection.size());

// Compute VH/VV ratio per image (speckle-filtered bands: VH_1, VV_1)
// Ratio increases with snow depth due to enhanced volume scattering
function computeRatio(image) {
  var vh  = image.select('VH_1').clip(AOI).clamp(-20, -7);
  var vv  = image.select('VV_1').clip(AOI).clamp(-20, -7);
  return vh.divide(vv).rename('ratio')
    .set('system:time_start', image.get('system:time_start'));
}

var ratioCollection = s1Collection.map(computeRatio);
var ratioList       = ratioCollection.toList(ratioCollection.size());

print('SAR ratio images available:', ratioList.size());

// Extract five acquisitions spanning the snow season
var r1 = ee.Image(ratioList.get(0));
var r2 = ee.Image(ratioList.get(1));
var r3 = ee.Image(ratioList.get(2));
var r4 = ee.Image(ratioList.get(3));
var r5 = ee.Image(ratioList.get(4));

// Cumulative sum of VH/VV ratios across all five acquisitions.
// Integrates the seasonal backscatter signal rather than relying
// on a single date — improves robustness against wet-snow events
// and temporal noise on individual scenes.
var siCumulative = r1.add(r2).add(r3).add(r4).add(r5).clamp(0.1, 1.9);
Map.addLayer(siCumulative, {}, 'SAR Cumulative SI (VH/VV ×5)');


// ── Snow Depth Model ─────────────────────────────────────────
// Semi-empirical model:
//
//   SD (cm) = [100 × SI_cumulative × a] / [1 − b × (NDVI − MSAVI)]
//
//   a = SD_COEFF_A (0.82) — empirical scaling coefficient
//   b = SD_COEFF_B (0.6)  — vegetation dampening factor
//
// The vegetation sub-index (NDVI − MSAVI) in the denominator
// corrects for canopy effects: higher residual vegetation signal
// reduces the denominator, boosting SD estimates proportionally.
// Model is restricted to optical snow mask (NDSI > threshold).

var sdNumerator   = siCumulative.multiply(100).multiply(SD_COEFF_A);
var sdDenominator = ee.Image(1).subtract(ee.Image(SD_COEFF_B).multiply(subIndex));
var sd            = sdNumerator.divide(sdDenominator).clip(AOI);

// NOTE — Alternative formulation with fixed vegetation fraction (FC):
// Uses a = 1.1 and hardcoded FC = 0.2 instead of the dynamic sub-index.
// Useful as a sensitivity check or where vegetation data is unavailable.
// Uncomment to compare:
//
// var sd_fixedFC = siCumulative.multiply(100).multiply(1.1)
//                   .divide(ee.Image(1).subtract(ee.Image(0.6).multiply(0.2)));
// Map.addLayer(sd_fixedFC.clamp(20, 350).clip(snowPoly), {palette: SD_PALETTE}, 'SD Fixed FC');


// ── Snow Depth Statistics ────────────────────────────────────
// Use 10th–90th percentile for robust visualisation stretch,
// avoiding influence of extreme outliers at snow margins.
print('── Computing snow depth statistics ──');

var sdMin = ee.Number(sd.reduceRegion({
  reducer: ee.Reducer.percentile([10]), geometry: snowPoly,
  scale: 100, maxPixels: 1e13, bestEffort: true
}).values().get(0)).round();

var sdMax = ee.Number(sd.reduceRegion({
  reducer: ee.Reducer.percentile([90]), geometry: snowPoly,
  scale: 100, maxPixels: 1e13, bestEffort: true
}).values().get(0)).round();

print('Snow Depth — 10th percentile (cm):', sdMin);
print('Snow Depth — 90th percentile (cm):', sdMax);

var sdMean = ee.Number(sd.reduceRegion({
  reducer: ee.Reducer.mean(), geometry: snowPoly,
  scale: 100, maxPixels: 1e13, bestEffort: true
}).values().get(0)).round();
print('Snow Depth — Mean (cm):', sdMean);


// ── Display Snow Depth ───────────────────────────────────────
Map.addLayer(
  sd.clip(snowPoly),
  { min: sdMin.getInfo(), max: sdMax.getInfo(), palette: SD_PALETTE },
  'Snow Depth (cm)'
);


// ── Legend ───────────────────────────────────────────────────
var colorBar = ui.Thumbnail({
  image : ee.Image.pixelLonLat().select(0),
  params: { bbox: [0, 0, 1, 0.1], dimensions: '100x10', format: 'png',
            min: 0, max: 1, palette: SD_PALETTE },
  style : { stretch: 'horizontal', margin: '0 8px', maxHeight: '24px' }
});

var legendLabels = ui.Panel({
  widgets: [
    ui.Label(sdMin,  { margin: '4px 8px' }),
    ui.Label('cm',   { margin: '4px 8px', textAlign: 'center', stretch: 'horizontal' }),
    ui.Label(sdMax,  { margin: '4px 8px' })
  ],
  layout: ui.Panel.Layout.flow('horizontal')
});

Map.add(ui.Panel([
  ui.Label({ value: 'Snow Depth (cm)', style: { fontWeight: 'bold', textAlign: 'center' } }),
  colorBar,
  legendLabels
]));


// ── Export SAR Ratio Images ───────────────────────────────────
// Exports each individual VH/VV ratio image for offline analysis.
[r1, r2, r3, r4, r5].forEach(function(img, i) {
  Export.image.toDrive({
    image         : img,
    description   : 'SAR_VH_VV_ratio_t' + (i + 1),
    fileNamePrefix: 'SAR_ratio_t' + (i + 1) + '_' + S1_START.slice(0, 7),
    folder        : EXPORT_FOLDER,
    scale         : EXPORT_SCALE,
    region        : AOI,
    crs           : 'EPSG:4326',
    maxPixels     : 1e13,
    fileFormat    : 'GeoTIFF'
  });
});


// ── Export Snow Depth ─────────────────────────────────────────
// Exported only within the optical snow mask (NDSI > threshold).
// Model is not physically valid outside snow-covered areas.
Export.image.toDrive({
  image         : sd.clip(snowPoly),
  description   : EXPORT_DESCRIPTION,
  fileNamePrefix: EXPORT_FILENAME,
  folder        : EXPORT_FOLDER,
  scale         : EXPORT_SCALE,
  region        : AOI,
  maxPixels     : 1e13,
  fileFormat    : 'GeoTIFF',
  formatOptions : { cloudOptimized: true }
});

print('── Analysis complete — check Tasks tab to export ──');
