Application.prototype.renderHypothesisListAt = function(container) {
  var chunks, hypothesis, row, template;
  chunks = [];
  template = (function() {
    var i, len, ref, results;
    ref = _.sortBy(_.values(this.hypotheses), function(o) {
      return o.name;
    });
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      hypothesis = ref[i];
      row = $("<div data-id=" + hypothesis.id + " class=''> <span class='handle'>" + (hypothesis.name.replace('_', ' ')) + "</span> <small class='right'>" + (hypothesis.hypothesis_type != null ? hypothesis.hypothesis_type : '-') + " <a class='remove-hypothesis mdi-content-clear'></a></small></div>");
      chunks.push(row);
      results.push(row.find(".remove-hypothesis").click(function(event) {
        var handle;
        return handle = $(this);
      }));
    }
    return results;
  }).call(this);
  return $(container).html(chunks);
};

Application.initializers.push(function() {
  return this.on("render-hypotheses", (function(_this) {
    return function() {
      return _this.renderHypothesisListAt(".hypothesis-list");
    };
  })(this));
});

//# sourceMappingURL=hypothesis-ui.js.map
