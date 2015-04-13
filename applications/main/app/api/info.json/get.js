module.exports = function(client, callback) {
  'use strict';

  var http = require('http');
  var zlib = require('zlib');

  var STATION_LATITUDE = 60;
  var STATION_LONGITUDE = 70;

  var query = '?latitude=' + STATION_LATITUDE +
              '&longitude=' + STATION_LONGITUDE;

  http.get('http://192.168.0.127/info' + query, function(res) {
    var chunks = [];

    res.on('data', function(chunk) {
      chunks.push(chunk);
    });

    res.on('end', function() {
      var buffer = Buffer.concat(chunks);
      
      if (res.headers['content-encoding'] == 'gzip') {
        processGzip(buffer);
      } else {
        process(buffer);
      }

    });

  }).on('error', function(e) {

    client.context.data = {
      error: e.message,
      explanation: ""
    };
    
    if (e.message == "connect ETIMEDOUT" ||
        e.message == "connect EHOSTUNREACH" ||
        e.message == "connect ECONNREFUSED") {
        
      client.context.data.explanation =
      "Targets generation server is not responding";
    }
    
    callback();

  });

  function processGzip(buffer) {
    zlib.gunzip(buffer, function(err, decoded) {

      var coords = JSON.parse(decoded.toString());
      // need processing
      client.context.data = {
        latitude: STATION_LATITUDE,
        longitude: STATION_LONGITUDE,
        coords: coords
      };
      callback();

    });
  }
  
  function process(buffer) {
    var coords = JSON.parse(buffer.toString());
    // need processing
    client.context.data = {
      latitude: STATION_LATITUDE,
      longitude: STATION_LONGITUDE,
      coords: coords
    };
    callback();
  }
};
