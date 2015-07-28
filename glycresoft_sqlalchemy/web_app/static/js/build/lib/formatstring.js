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
      var err, v;
      if (name === '') {
        name = keys[i];
        i++;
      }
      try {
        v = JSON.stringify(data[name]);
        if (v.length > 1) {
          v = v.slice(1, -1);
        }
        return v;
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
