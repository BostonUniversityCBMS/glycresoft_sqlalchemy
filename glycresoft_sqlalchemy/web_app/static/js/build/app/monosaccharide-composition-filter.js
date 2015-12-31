var MonosaccharideFilter;

MonosaccharideFilter = (function() {
  function MonosaccharideFilter(parent, residueNames, rules) {
    if (rules == null) {
      if (GlycReSoft.settings.monosaccharide_filters == null) {
        GlycReSoft.settings.monosaccharide_filters = {};
      }
      rules = GlycReSoft.settings.monosaccharide_filters;
    }
    this.container = $("<div></div>").addClass("row");
    $(parent).append(this.container);
    this.residueNames = residueNames;
    this.rules = rules;
  }

  MonosaccharideFilter.prototype.makeFilterWidget = function(residue) {
    var rendered, rule, sanitizeName, self, template;
    rule = this.rules[residue];
    if (rule == null) {
      rule = {
        minimum: 0,
        maximum: 10,
        include: true
      };
      this.rules[residue] = rule;
    }
    residue.name = residue;
    residue.sanitizeName = sanitizeName = residue.replace(/[\(\),]/g, "_");
    template = "<span class=\"col s2 monosaccharide-filter\" data-name='" + residue + "'>\n    <p style='margin: 0px; margin-bottom: -10px;'>\n        <input type=\"checkbox\" id=\"" + sanitizeName + "_include\" name=\"" + sanitizeName + "_include\"/>\n        <label for=\"" + sanitizeName + "_include\"><b>" + residue + "</b></label>\n    </p>\n    <p style='margin-top: 0px; margin-bottom: 0px;'>\n        <input id=\"" + sanitizeName + "_min\" type=\"number\" placeholder=\"Minimum " + residue + "\" style='width: 45px;' min=\"0\"\n               value=\"" + rule.minimum + "\" max=\"" + rule.maximum + "\" name=\"" + sanitizeName + "_min\"/> : \n        <input id=\"" + sanitizeName + "_max\" type=\"number\" placeholder=\"Maximum " + residue + "\" style='width: 45px;' min=\"0\"\n               value=\"" + rule.maximum + "\" name=\"" + sanitizeName + "_max\"/>\n    </p>\n</span>";
    self = this;
    rendered = $(template);
    rendered.find("#" + sanitizeName + "_min").change(function() {
      rule.minimum = parseInt($(this).val());
      return self.changed();
    });
    rendered.find("#" + sanitizeName + "_max").change(function() {
      rule.maximum = parseInt($(this).val());
      return self.changed();
    });
    rendered.find("#" + sanitizeName + "_include").prop("checked", rule.include).click(function() {
      rule.include = $(this).prop("checked");
      return self.changed();
    });
    return rendered;
  };

  MonosaccharideFilter.prototype.render = function() {
    var i, len, ref, residue, results, widget;
    ref = this.residueNames;
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      residue = ref[i];
      widget = this.makeFilterWidget(residue);
      results.push(this.container.append(widget));
    }
    return results;
  };

  MonosaccharideFilter.prototype.changed = _.debounce((function() {
    return GlycReSoft.emit("update_settings");
  }), 1000);

  return MonosaccharideFilter;

})();

//# sourceMappingURL=monosaccharide-composition-filter.js.map
