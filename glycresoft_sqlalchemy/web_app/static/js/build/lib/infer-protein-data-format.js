var getProteinName, getProteinNamesFromMzIdentML, identifyProteomicsFormat;

identifyProteomicsFormat = function(file, callback) {
  var isMzidentML, reader;
  isMzidentML = function(lines) {
    var i, len, line;
    for (i = 0, len = lines.length; i < len; i++) {
      line = lines[i];
      if (/mzIdentML/.test(line)) {
        return true;
      }
    }
    return false;
  };
  reader = new FileReader();
  reader.onload = function() {
    var lines, proteomicsFileType;
    lines = this.result.split("\n");
    console.log(lines);
    proteomicsFileType = "fasta";
    if (isMzidentML(lines)) {
      proteomicsFileType = "mzIdentML";
    }
    return callback(proteomicsFileType);
  };
  return reader.readAsText(file.slice(0, 100));
};

getProteinName = function(line) {
  return line.split("_", 2)[1];
};

getProteinNamesFromMzIdentML = function(file, callback) {
  var chunksize, fr, offset, proteins, seek;
  fr = new FileReader();
  chunksize = 1024 * 8;
  offset = 0;
  proteins = {};
  fr.onload = function() {
    var i, len, line, lines, name;
    lines = this.result.split("\n");
    for (i = 0, len = lines.length; i < len; i++) {
      line = lines[i];
      if (/<ProteinDetectionHypothesis/i.test(line)) {
        name = getProteinName(line);
        console.log(name);
        proteins[name] = true;
      }
    }
    return seek();
  };
  fr.onerror = function(error) {
    return console.log(error);
  };
  seek = function() {
    if (offset >= file.size) {
      return callback(Object.keys(proteins));
    } else {
      fr.readAsText(file.slice(offset, offset + chunksize));
      return offset += chunksize / 2;
    }
  };
  return seek();
};

//# sourceMappingURL=infer-protein-data-format.js.map
