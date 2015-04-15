module.exports = function(client, callback) {
  'use strict';

  var http = require('http');
  var zlib = require('zlib');

  var STATION_LATITUDE = 44.93;
  var STATION_LONGITUDE = 40.98;

  var query = '?latitude=' + STATION_LATITUDE +
              '&longitude=' + STATION_LONGITUDE;
              
  console.log('...request');

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
      process(decoded);
    });
  }
  
  function process(buffer) {
    var coords = JSON.parse(buffer.toString());
    // need processing
    
    coords.latitude = STATION_LATITUDE;
    coords.longitude = STATION_LONGITUDE;
    
    client.context.data = coords;
    callback();
  }
};
