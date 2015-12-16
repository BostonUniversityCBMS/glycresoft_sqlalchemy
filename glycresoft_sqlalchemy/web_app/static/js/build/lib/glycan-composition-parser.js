var GlycanComposition;

GlycanComposition = (function() {
  function GlycanComposition(string) {
    this.__order = [];
    this.map = {};
    this.parse(string);
  }

  GlycanComposition.prototype.parse = function(string) {
    var i, len, name, number, part, parts, ref, results;
    parts = string.slice(1, -1).split("; ");
    results = [];
    for (i = 0, len = parts.length; i < len; i++) {
      part = parts[i];
      ref = part.split(":"), name = ref[0], number = ref[1];
      this.__order.push(name);
      results.push(this.map[name] = parseInt(number));
    }
    return results;
  };

  GlycanComposition.prototype.format = function(colorSource) {
    var color, name, number, parts, ref, template;
    parts = [];
    ref = this.map;
    for (name in ref) {
      number = ref[name];
      console.log(name, number);
      if (name === '__order') {
        continue;
      }
      color = colorSource.get(name);
      template = "<span class='monosaccharide-name' style='background-color:" + color + "; padding: 2px;border-radius:2px;'>" + name + " " + number + "</span>";
      parts.push(template);
    }
    return parts.join(' ');
  };

  return GlycanComposition;

})();

//# sourceMappingURL=glycan-composition-parser.js.map
