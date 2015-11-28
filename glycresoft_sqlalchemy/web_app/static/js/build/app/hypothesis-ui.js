Application.prototype.renderHypothesisListAt = function(container) {
  var chunks, hypothesis, i, len, ref, row, self, template;
  chunks = [];
  template = '';
  self = this;
  ref = _.sortBy(_.values(this.hypotheses), function(o) {
    return o.name;
  });
  for (i = 0, len = ref.length; i < len; i++) {
    hypothesis = ref[i];
    row = $("<div data-id=" + hypothesis.id + " class=''> <span class='handle'>" + (hypothesis.name.replace('_', ' ')) + "</span> <small class='right'>" + (hypothesis.hypothesis_type != null ? hypothesis.hypothesis_type : '-') + " <a class='remove-hypothesis mdi-content-clear'></a></small></div>");
    chunks.push(row);
    row.click(function(event) {
      var handle, hypothesisId, layer;
      handle = $(this);
      hypothesisId = handle.attr("data-id");
      self.addLayer(ActionBook.viewHypothesis, {
        "hypothesis_id": hypothesisId
      });
      layer = self.lastAdded;
      return self.setShowingLayer(layer);
    });
    row.find(".remove-hypothesis").click(function(event) {
      var handle;
      return handle = $(this);
    });
  }
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
