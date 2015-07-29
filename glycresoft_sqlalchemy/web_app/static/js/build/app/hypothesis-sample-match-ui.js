Application.prototype.renderHypothesisSampleMatchListAt = function(container) {
  var chunks, hsm, row, self, template;
  chunks = [];
  template = (function() {
    var i, len, ref, results;
    ref = _.sortBy(_.values(this.hypothesisSampleMatches), function(o) {
      return o.name;
    });
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      hsm = ref[i];
      row = $("<div data-id=" + hsm.id + " class=''> <span class='handle'>" + (hsm.name.replace('_', ' ')) + "</span> <small class='right'>" + " <a class='remove-hsm mdi-content-clear'></a></small></div>");
      chunks.push(row);
      self = this;
      row.click(function(event) {
        var handle, id;
        handle = $(this);
        id = handle.attr('data-id');
        self.addLayer(ActionBook.viewDatabaseSearchResults, [id]);
        console.log(self.layers);
        console.log(self.lastAdded);
        return self.setShowingLayer(self.lastAdded);
      });
      results.push(row.find(".remove-hsm").click(function(event) {
        var handle;
        return handle = $(this);
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
