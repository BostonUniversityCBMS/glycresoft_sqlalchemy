Application.prototype.renderHypothesisSampleMatchListAt = function(container) {
  var chunks, hsm, row, self, template;
  chunks = [];
  template = (function() {
    var i, len, ref, results;
    ref = _.sortBy(_.values(this.hypothesisSampleMatches), function(o) {
      return o.id;
    });
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      hsm = ref[i];
      hsm.name = hsm.name != null ? hsm.name : "HypothesisSampleMatch:" + hsm.target_hypothesis.name + "@" + hsm.sample_run_name;
      row = $("<div data-id=" + hsm.id + " class='list-item clearfix'> <span class='handle'>" + hsm.id + ". " + (hsm.name.replace('_', ' ')) + "</span> <small class='right' style='display:inherit'> " + (hsm.hypothesis_sample_match_type.replace('HypothesisSampleMatch', '')) + " <a class='remove-hsm mdi-content-clear'></a> </small> </div>");
      chunks.push(row);
      self = this;
      row.click(function(event) {
        var handle, id;
        handle = $(this);
        id = handle.attr('data-id');
        self.addLayer(ActionBook.viewDatabaseSearchResults, {
          hypothesis_sample_match_id: id
        });
        console.log(self.layers);
        console.log(self.lastAdded);
        self.context["hypothesis_sample_match_id"] = id;
        return self.setShowingLayer(self.lastAdded);
      });
      results.push(row.find(".remove-hsm").click(function(event) {
        var handle;
        handle = $(this);
        return console.log("Removal of HypothesisSampleMatch is not implemented.");
      }));
    }
    return results;
  }).call(this);
  return $(container).html(chunks);
};

Application.initializers.push(function() {
  return this.on("render-hypothesis-sample-matches", (function(_this) {
    return function() {
      return _this.renderHypothesisSampleMatchListAt(".hypothesis-sample-match-list");
    };
  })(this));
});

//# sourceMappingURL=hypothesis-sample-match-ui.js.map
