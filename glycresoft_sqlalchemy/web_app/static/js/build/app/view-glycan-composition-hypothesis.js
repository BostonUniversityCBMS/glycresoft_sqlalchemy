var viewGlycanCompositionHypothesis;

viewGlycanCompositionHypothesis = function(hypothesisId) {
  var currentPage, detailModal, displayTable, setup, setupGlycanCompositionTablePageHandler, updateCompositionTablePage;
  detailModal = void 0;
  displayTable = void 0;
  currentPage = 1;
  setup = function() {
    displayTable = $("#composition-table-container");
    return updateCompositionTablePage(1);
  };
  setupGlycanCompositionTablePageHandler = function(page) {
    if (page == null) {
      page = 1;
    }
    $('.display-table tbody tr').click(function() {});
    $(':not(.disabled) .next-page').click(function() {
      return updateCompositionTablePage(page + 1);
    });
    $(':not(.disabled) .previous-page').click(function() {
      return updateCompositionTablePage(page - 1);
    });
    return $('.pagination li :not(.active)').click(function() {
      var nextPage;
      nextPage = $(this).attr("data-index");
      if (nextPage != null) {
        nextPage = parseInt(nextPage);
        return updateCompositionTablePage(nextPage);
      }
    });
  };
  updateCompositionTablePage = function(page) {
    var url;
    if (page == null) {
      page = 1;
    }
    url = "/view_glycan_composition_hypothesis/" + hypothesisId + "/" + page;
    console.log(url);
    return GlycReSoft.ajaxWithContext(url).success(function(doc) {
      currentPage = page;
      displayTable.html(doc);
      return setupGlycanCompositionTablePageHandler(page);
    });
  };
  return setup();
};

//# sourceMappingURL=view-glycan-composition-hypothesis.js.map
