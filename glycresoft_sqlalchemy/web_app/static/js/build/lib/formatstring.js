(function() {

  /*
  Implements {named} replacements similar to the simple format() method of strings from Python
   */
  String.prototype.format = function() {
    var data, i, keys, res;
    data = arguments;
    i = 0;
    keys = Object.keys(arguments);
    if (arguments.length === 1 && typeof arguments[0] === 'object') {
      data = arguments[0];
      keys = Object.keys(arguments);
    }
    res = this.replace(/\{([^\}]*)\}/g, function(placeholder, name, position) {
      var err;
      if (name === '') {
        name = keys[i];
        i++;
      }
      try {
        return JSON.stringify(data[name]).slice(1, -1);
      } catch (_error) {
        err = _error;
        console.log(err, name, data);
        return void 0;
      }
    });
    return res;
  };
})();

//# sourceMappingURL=formatstring.js.map
