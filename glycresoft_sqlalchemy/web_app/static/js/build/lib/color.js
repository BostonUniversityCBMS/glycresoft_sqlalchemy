var ColorManager;

ColorManager = (function() {
  function ColorManager(obj) {
    this.store = obj;
  }

  ColorManager.prototype.get = function(name) {
    var b, g, r, rgb, string;
    try {
      rgb = this.store[name];
      r = rgb[0], g = rgb[1], b = rgb[2];
      string = "rgba(" + (r * 255) + "," + (g * 255) + "," + (b * 255) + ",0.5)";
      return string;
    } catch (_error) {
      return this.update;
    }
  };

  ColorManager.prototype.update = function(callback) {
    if (callback == null) {
      callback = function() {
        return {};
      };
    }
    return $.getJSON("/api/colors").success((function(_this) {
      return function(data) {
        _this.store = data;
        return callback(_this);
      };
    })(this));
  };

  return ColorManager;

})();

//# sourceMappingURL=color.js.map
