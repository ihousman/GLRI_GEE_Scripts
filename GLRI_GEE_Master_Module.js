//L5,L7,L8  simplelandsatcomposite-6/19/14

//Code written by:
//Cloud/Shadow masking- Carson Stam
//Image normalization- Kim McCallum
//Change detection and module integration- Ian Housman

//RedCastle Resources, Inc.

//Working onsite at:
//USDA Forest Service
//Remote Sensing Applications Center (RSAC)
//2222 West 2300 South
//Salt Lake City, UT 84119
//Office: (801) 975-3366
//Email: ihousman@fs.fed.us
//RSAC FS Intranet website: http://fsweb.rsac.fs.fed.us/
//RSAC FS Internet website: http://www.fs.fed.us/eng/rsac/
//////////////////////////////////////////////////////////////////
var range = function(start, stop, step){
    if (typeof stop=='undefined'){
        // one param defined
        stop = start;
        start = 0;
    };
    if (typeof step=='undefined'){
        step = 1;
    };
    if ((step>0 && start>=stop) || (step<0 && start<=stop)){
        return [];
    };
    var result = [];
    for (var i=start; step>0 ? i<stop : i>stop; i+=step){
        result.push(i);
    };
    return result;
};
var stringify = function(in_no){
  return in_no.toString();
}
///////////////////////////////////////////////////////////////
//Specify the dates
  var collection_dict = {L8: 'LC8_L1T_TOA',
                         L7: 'L7_L1T_TOA',
                         L5: 'L5_L1T_TOA'};
                       
  var sensor_year_dict = {L8: range(2013,2020),
                         L7: range(2009,2013),
                         L5: range(1985, 2009)};
                         
  var shadowPop_sensor_year_dict = {L8: range(2013,2020),
                         L7: range(1999,2020),
                         L5: range(1985, 2008)};
                         
                         
  var sensor_band_dict ={L8 : [1,2,3,4,5,9,6],
                        L7 : [0,1,2,3,4,5,7],
                        L5 : [0,1,2,3,4,5,6]};
                        
var sensorShadowSumDict ={L8 : [4,5,6],
                          L7 : [3,4,7],
                          L5 : [3,4,6]};
var sensorIncompleteSumDict ={L8 : [1,2,3,4,5,6,7],
                          L7 : [0,1,2,3,4,5,7],
                          L5 : [0,1,2,3,4,5,6]};
var sensorDarkThresholdsDict ={L8 : [0.3,0.1,0.18],
                              L7 : [0.15,0.12,0.1],
                              L5 : [0.15,0.12,0.1]};
                              
var possible_sensors = ['L8','L7','L5'];
var spacecraft_dict = {'Landsat5': 'L5','Landsat7': 'L7', 'LANDSAT_8': 'L8'}

var STD_NAMES = ['blue', 'green', 'red', 'nir', 'swir1', 'temp','swir2'];

var years = [1990];//,1990, 1995,2000,2005,2010];//,1995,2000,2005,2010,2013];
var start_julian = 190;
var end_julian = 258;

var composite_year_no = 2;
var percentile =50;
var cloud_thresh = 10;
var shadowThresh = -1;
var shadowPopulationYearLookBack = 3;
var shadowPopulationYearLookForward = 5;

//Buffering options for clouds masking
//Likely only needed for single-image masking
var bufferCS = false;
var bufferCS_contract = 1;
var bufferCS_expand = 3;

//Decides whether shadow pop needs to be computed across a variable names years
var run_multi_date = true;

var flat_threshold = 0;
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//Set up flat areas for water masking
var flat = ee.call('Terrain', ee.Image('srtm90_v4')).select('slope').lte(flat_threshold);
///////////////////////////////////////////////////////

//var all_bands = '';
var all_bands_list = [];
var i = 0;
///////////////////////////////////////////////////////
//Funds the sensor for a given year
//Only year is a required parameter

var sensor_year_finder = function(year, syDict,return_list){
  year = parseInt(year);
  var sensor = 'L5';
  var sensor_list = [];
  
  if(typeof syDict === 'undefined'){syDict = sensor_year_dict};
  
  if(typeof return_list === 'undefined'){return_list = false};

  for (var i in possible_sensors){
    var poss_sensor = possible_sensors[i];
    if (syDict[poss_sensor].indexOf(year)!= -1){
      sensor= poss_sensor;
      sensor_list.push(poss_sensor);
     
      }
  }
    if(return_list === false){return sensor}
    else{return sensor_list};
    
};

//////////////////////////////////////////////////////////////////////////
//Functions for shadow masking
//Author: Carson Stam
//Adapted by: Ian Housman
//////////////////////////////////////////////////////////////////////////

var sum_fun = function(image,sensor){
  return image.select(sensorShadowSumDict[sensor]).reduce(ee.Reducer.sum());
}
//////////////////////////////////////////////////////////////////////////
//Function for masking out incomplete landsat
//Author: Carson Stam
//Adapted by: Ian Housman
//////////////////////////////////////////////////////////////////////////

var maskIncomplete = function(image, sensorName)
  {
    var band_list = sensorIncompleteSumDict[sensorName];
    var imageWhere = image.where(image.select([band_list[0]]).gt(-0.5).or(image.select([band_list[1]]).gt(-0.5)).or(image.select([band_list[2]]).gt(-0.5)).or(image.select([band_list[3]]).gt(-0.5)).or(image.select([band_list[4]]).gt(-0.5)).or(image.select([band_list[5]]).gt(-0.5)).or(image.select([band_list[6]]).gt(-0.5)),10);
    return image.mask(imageWhere.select([1]).eq(10));
  };

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Function for masking out dark pixels
//Author: Carson Stam
//Adapted by: Ian Housman
//////////////////////////////////////////////////////////////////////////
var maskDark = function(image,sensor){
  var SumBands = sensorShadowSumDict[sensor];
  var SumThreshs = sensorDarkThresholdsDict[sensor];
  var darkOut = image.where
    (
    image.select(SumBands[2]).gte(SumThreshs[2])
    .or(image.select(SumBands[0]).gte(SumThreshs[0]))
    .or(image.select(SumBands[1]).gte(SumThreshs[1])),3
    );
  return darkOut.select([1]);
  };


/////////////////////////////////////////////////////////////
//Shadow population getter
//Author: Carson Stam
//Adapted by: Ian Housman
/////////////////////////////////////////////////////
var shadowPopulation = function(year)
  {
    print('Shadow popping', year);
    //Get bounding year window
    //if(typeof lookBackPeriod == 'undefined'){var lookBackPeriod = 3};
    var y1 = year - shadowPopulationYearLookBack;
    var y2 = year + shadowPopulationYearLookForward;
    
    //Find sensors bounding year window intersects
    var sD = new Date('1/1/' + y1.toString());
    var eD = new Date('12/31/' + y2.toString());
    
    var sensors = sensor_year_finder(year,shadowPop_sensor_year_dict,true);
    var sensor = "L5";
    
    if(sensors.indexOf(sensor) != -1)
      {
      //Get the collection
      var temporalStackL5 = ee.ImageCollection(collection_dict[sensor])
                      .filterDate(sD, eD)
                      .filter(ee.Filter.calendarRange(start_julian, end_julian));
      //Mask out and sum up the collection
      temporalStackL5 = temporalStackL5.map(function(maskOUT){return maskIncomplete(maskOUT,sensor)});
      temporalStackL5 = temporalStackL5.map(function(tsl5){return sum_fun(tsl5,sensor)});
      var outStack = temporalStackL5;
      }
    sensor = "L7";
    if(sensors.indexOf(sensor) != -1)
      {
      //Get the collection
      var temporalStackL7 = ee.ImageCollection(collection_dict[sensor])
                      .filterDate(sD, eD)
                      .filter(ee.Filter.calendarRange(start_julian, end_julian));
      //Mask out and sum up the collection
      temporalStackL7 = temporalStackL7.map(function(maskOUT){
                                        var mOut = maskIncomplete(maskOUT,sensor);
                                        return sum_fun(mOut,sensor)});
      if(outStack === undefined){var outStack = temporalStackL7}
      else{outStack = ee.ImageCollection(outStack.merge(temporalStackL7))};
      }
    sensor = "L8";
    if(sensors.indexOf(sensor) != -1)
      {
      //Get the collection
      var temporalStackL8 = ee.ImageCollection(collection_dict[sensor])
                      .filterDate(sD, eD)
                      .filter(ee.Filter.calendarRange(start_julian, end_julian));
      //Mask out and sum up the collection
      temporalStackL8 = temporalStackL8.map(function(maskOUT){
                                        var mOut = maskIncomplete(maskOUT,sensor);
                                        return sum_fun(mOut,sensor)});
      if(outStack === undefined){var outStack = temporalStackL8}
      else{outStack = ee.ImageCollection(outStack.merge(temporalStackL8))};                                  
      }
  
    return outStack;
  };
//////////////////////////////////////////////////////////////////////////

//Create array of shadow Populations
//////////////////////////////////////////////////////////////////////////
if(run_multi_date === true){
var shadowYear1 = years.sort()[0]-shadowPopulationYearLookBack;
var shadowYear2 = years.sort()[years.length-1];
shadowYear2 = shadowYear2 + shadowPopulationYearLookForward;
var yearArray = range(shadowYear1,shadowYear2);
//Create array of shadowpopulations
var shadowPopArray = yearArray.map(shadowPopulation);
}

//////////////////////////////////////////////////////////////////////////
//Functions for cloud masking
//////////////////////////////////////////////////////////////////////////
var cloudScore = function(cloudImage){
  cloudImage = ee.Algorithms.Landsat.simpleCloudScore(cloudImage);
  return cloudImage.select('cloud');
};
//////////////////////////////////////////////////////////////////////////
//Bust out l8 clouds using cloud score instead of BQA
var bust_clouds = function(image) {
  
  image = ee.Algorithms.Landsat.simpleCloudScore(image);
  var quality = image.select('cloud');
  var cloud01 = quality.gt(cloud_thresh);
  var maskedImage = image.mask().and(cloud01.not());
  image = image.mask(maskedImage)
  return image;
};
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//Function for applying zcore and masking cloud shadow
//Author: Carson Stam
//Adapted by: Ian Housman
//////////////////////////////////////////////////////////////////////////
var shadowMasking = function(image, year, yearArray_local,shadowPopArray_local)
  {
    
    if(typeof yearArray_local != undefined){var yearArray = yearArray_local};
    if(typeof shadowPopArray_local != undefined){var shadowPopArray = shadowPopArray_local};
    
    print('yeararray',yearArray);
    print('shadowpoparray',shadowPopArray);
    var sensor = sensor_year_finder(year);
    var yearIndex = yearArray.indexOf(year);
    
    var landsatShadowPop = ee.ImageCollection(shadowPopArray[yearIndex]);
    //var mean = landsatShadowPop.reduce(ee.Reducer.mean());
    var mean = landsatShadowPop.mean();
    var sigma = landsatShadowPop.reduce(ee.Reducer.stdDev()); 

    var darkAreas = maskDark(image,sensor);
    var ImageSum = sum_fun(image,sensor);
    var zscore = ImageSum.subtract(mean).divide(sigma);
    
    if(bufferCS)
      {
        var imageCloudScore = cloudScore(image);
        var cloudShadowMask = (zscore.lt(shadowThresh).and(darkAreas.neq(3)).or(imageCloudScore.gt(cloud_thresh)));
        var cloudShadowMaskFilter = cloudShadowMask.focal_min(bufferCS_contract).focal_max(bufferCS_expand);
        var outImage =  image.mask(cloudShadowMaskFilter.eq(0));
      }
    else
      {
        image = bust_clouds(image);
        var imageWhere = image.where((zscore.lt(shadowThresh).and(darkAreas.neq(3))),5);
        var outImage =  image.mask(imageWhere.neq(5));
      }
    return maskIncomplete(outImage,sensor);
    
  };

var image_getter = function(year){
  //Find collection name from year
  var sensor = sensor_year_finder(year);
  var collection_name = collection_dict[sensor];
  var bands = sensor_band_dict[sensor];

  //Set up dates
  var startYear = year.toString();
  var endYear = (year + composite_year_no).toString();
  var startDate = new Date('1/1/' + startYear);
  var endDate = new Date('12/31/' + endYear);
  print([startYear,endYear, collection_name,bands]);
  //print(i);
  
    
  var collection =  ee.ImageCollection(collection_name)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.calendarRange(start_julian, end_julian));
    //.filterMetadata('WRS_PATH', 'EQUALS', 41)
    //.filterMetadata('WRS_ROW', 'EQUALS', 26);
  
  var set_date_time_add_indices = function(in_image){
    in_image = in_image.addBands(in_image.normalizedDifference(['nir', 'red']).select([0],['ndvi']));
  var time_start = new Date(year,1,1);
  var time_end = new Date(year+composite_year_no,12,31);
  
  //Set up a water mask
  var mndwi =in_image.normalizedDifference(['green', 'swir1']).select([0],['mndwi']);
  var water_mask = mndwi.mask(mndwi.gt(0.3).and(flat.eq(1))).select([0],['water_mask']);
  in_image = in_image.addBands(water_mask);
  
  //Set some properties for book keeping down the road
  in_image = in_image.set({
      'system:time_start': time_start.valueOf(),
      'system:time_end': time_end.valueOf(),
      'system:year_start': year,
      'system:year_end':year+composite_year_no
  });
  return in_image;
  }

   
  //Cloud and shadow masking
  //Shadow masking written by: Carson Stam
  //var shadowPop = shadowPopulation(year); 
  //addToMap(shadowPop,{},'ShadowPop_' + String(year));
  var shadow_masked = collection.map(function(ti){return shadowMasking(ti,year)}).reduce(ee.Reducer.percentile([percentile])).select(bands, STD_NAMES);
  shadow_masked = set_date_time_add_indices(shadow_masked);
  var busted = collection.map(bust_clouds).reduce(ee.Reducer.percentile([percentile])).select(bands, STD_NAMES);
  busted = set_date_time_add_indices(busted);

  all_bands_list.push(shadow_masked);
  addToMap(busted, {'max': 0.6, 'gamma': 0.8,  'bands':'swir1,nir,red'},'not_shadow_masked_'+startYear + '_' + endYear + '_' + collection_name, false);
  addToMap(shadow_masked, {'max': 0.6, 'gamma': 0.8,  'bands':'swir1,nir,red'},'shadow_masked_' +startYear + '_' + endYear + '_' + collection_name, true);
  addToMap(shadow_masked.select('water_mask'),{'palette': '0000FF,0000FF'},'water_shadow_masked_' +startYear + '_' + endYear + '_' + collection_name, false)
  i = i +1;
};

var singleImageMasking = function(img){
  
  var info = img.getInfo().properties;
  var year = parseInt(info.DATE_ACQUIRED.split("-",1)[0]);
  var spacecraft = info.SPACECRAFT_ID;
  var sensor = spacecraft_dict[spacecraft];
  var collection_name = collection_dict[sensor];
  var bands = sensor_band_dict[sensor];
  
  var shadowYear1 = year-shadowPopulationYearLookBack;
  var shadowYear2 = year + shadowPopulationYearLookForward;
  var yearArray = range(shadowYear1,shadowYear2);
 
//Create array of shadowpopulations
  var shadowPopArray = yearArray.map(shadowPopulation);
  
  print(sensor);
  var maskedImage = shadowMasking(img,year,yearArray,shadowPopArray).select(bands,STD_NAMES);
  print(maskedImage.getInfo());
  //addToMap(maskedImage.select(bands));

  
  //addToMap(img.select([bands[6],bands[3],bands[2]]));
  
}
//Single image example
//var img = ee.Image('L7_L1T_TOA/LE70410262012187EDC00');
//var img = ee.Image('LC8_L1T_TOA/LC80410262013197LGN00');
var img = ee.Image('L5_L1T_TOA/LT50410262009266PAC01')
//var sensor = 'L8';

//singleImageMasking(img);

//Map the image_getter function across all years to get 
//an image collection containing composites
years.map(image_getter);

//Cast the list to an image collection
//var all_bands = ee.ImageCollection(all_bands_list);//.select(range(0,years.length),years.map(stringify));


//var trends=all_bands.map(function (image) {
//  return image.select(['nir'])}).formaTrend().select('long-trend');//.mask(hansen_image_masked);//.clip(fc);

//addToMap(trends.mask(trends.lt(-0.001)), {'max': 0, 'min': -0.02, 'palette': 'FF0000,FFFF00'},'trend',false);



