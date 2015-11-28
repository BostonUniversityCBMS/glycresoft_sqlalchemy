var viewGlycopeptideCompositionHypothesis;

viewGlycopeptideCompositionHypothesis = function(hypothesisId) {
  var currentPage, displayTable, proteinContainer, proteinId, setup, setupGlycopeptideCompositionTablePageHandler, updateCompositionTablePage, updateProteinChoice;
  displayTable = void 0;
  currentPage = 1;
  proteinContainer = void 0;
  proteinId = void 0;
  setup = function() {
    proteinContainer = $("#protein-container");
    $('.protein-list-table tbody tr').click(updateProteinChoice);
    return updateProteinChoice.apply($('.protein-list-table tbody tr'));
  };
  setupGlycopeptideCompositionTablePageHandler = function(page) {
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
  updateProteinChoice = function() {
    var handle, id, url;
    handle = $(this);
    proteinId = id = handle.attr('data-target');
    proteinContainer.html("<div class=\"progress\"><div class=\"indeterminate\"></div></div>").fadeIn();
    url = "/view_glycopeptide_composition_hypothesis/protein_view/" + proteinId;
    return $.post(url, {
      "settings": GlycReSoft.settings,
      "context": GlycReSoft.context
    }).success(function(doc) {
      proteinContainer.hide();
      proteinContainer.html(doc).fadeIn();
      GlycReSoft.context["current_protein"] = id;
      displayTable = $("#display-table-container");
      return updateCompositionTablePage(1);
    }).error(function(error) {
      return console.log(arguments);
    });
  };
  updateCompositionTablePage = function(page) {
    var url;
    if (page == null) {
      page = 1;
    }
    url = "/view_glycopeptide_composition_hypothesis/protein_view/" + proteinId + "/" + page;
    console.log(url);
    return GlycReSoft.ajaxWithContext(url).success(function(doc) {
      currentPage = page;
      displayTable.html(doc);
      return setupGlycopeptideCompositionTablePageHandler(page);
    });
  };
  return setup();
};

//# sourceMappingURL=view-glycopeptide-composition-hypothesis.js.map
