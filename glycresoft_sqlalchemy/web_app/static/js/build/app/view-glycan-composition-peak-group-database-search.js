var doZoom, viewGlycanCompositionPeakGroupingDatabaseSearchResults;

doZoom = function(selector) {
  var svg, zoom;
  svg = d3.select(selector);
  zoom = function() {
    return svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
  };
  return d3.select(selector).call(d3.behavior.zoom().scaleExtent([1, 8]).on("zoom", zoom));
};

viewGlycanCompositionPeakGroupingDatabaseSearchResults = function() {
  var currentPage, downloadCSV, glycanDetailsModal, glycanTable, setup, setupGlycanCompositionTablePageHandlers, showGlycanCompositionDetailsModal, unload, updateGlycanCompositionTablePage, updateView;
  glycanDetailsModal = void 0;
  glycanTable = void 0;
  currentPage = 1;
  setup = function() {
    updateView();
    return $("#save-csv-file").click(downloadCSV);
  };
  setupGlycanCompositionTablePageHandlers = function(page) {
    if (page == null) {
      page = 1;
    }
    $('.glycan-match-row').click(showGlycanCompositionDetailsModal);
    $(':not(.disabled) .next-page').click(function() {
      return updateGlycanCompositionTablePage(page + 1);
    });
    $(':not(.disabled) .previous-page').click(function() {
      return updateGlycanCompositionTablePage(page - 1);
    });
    return $('.pagination li :not(.active)').click(function() {
      var nextPage;
      nextPage = $(this).attr("data-index");
      if (nextPage != null) {
        nextPage = parseInt(nextPage);
        return updateGlycanCompositionTablePage(nextPage);
      }
    });
  };
  updateGlycanCompositionTablePage = function(page) {
    var url;
    if (page == null) {
      page = 1;
    }
    url = "/view_database_search_results/glycan_composition_match_table/" + page;
    console.log(url);
    return GlycReSoft.ajaxWithContext(url).success(function(doc) {
      currentPage = page;
      glycanTable.html(doc);
      return setupGlycanCompositionTablePageHandlers(page);
    });
  };
  updateView = function() {
    var handle;
    handle = $(this);
    $("#content-container").html("<div class=\"progress\"><div class=\"indeterminate\"></div></div>").fadeIn();
    return GlycReSoft.ajaxWithContext('/view_database_search_results/results_view/').success(function(doc) {
      var tabs;
      $('#content-container').hide();
      $('#content-container').html(doc).fadeIn();
      tabs = $('ul.tabs');
      tabs.tabs();
      if (GlycReSoft.context['view-active-tab'] !== void 0) {
        console.log(GlycReSoft.context['view-active-tab']);
        $('ul.tabs').tabs('select_tab', GlycReSoft.context['view-active-tab']);
      } else {
        $('ul.tabs').tabs('select_tab', 'glycome-overview');
      }
      $('.indicator').addClass('indigo');
      $('ul.tabs .tab a').click(function() {
        return GlycReSoft.context['view-active-tab'] = $(this).attr('href').slice(1);
      });
      glycanDetailsModal = $('#glycan-detail-modal');
      glycanTable = $("#glycan-table");
      return updateGlycanCompositionTablePage(1);
    }).error(function(error) {
      return console.log(arguments);
    });
  };
  showGlycanCompositionDetailsModal = function() {
    var handle, id;
    handle = $(this);
    id = handle.attr('data-target');
    console.log(id);
    return PartialSource.glycanCompositionDetailsModal({
      "id": id
    }, function(doc) {
      glycanDetailsModal.find('.modal-content').html(doc);
      $(".lean-overlay").remove();
      return glycanDetailsModal.openModal();
    });
  };
  unload = function() {
    return GlycReSoft.removeCurrentLayer();
  };
  downloadCSV = function() {
    var handle, id;
    handle = $(this);
    id = handle.attr('data-target');
    return $.ajax("/view_database_search_results/export_csv/" + id, {
      data: JSON.stringify({
        "context": GlycReSoft.context,
        "settings": GlycReSoft.settings
      }),
      contentType: "application/json",
      type: 'POST'
    });
  };
  return setup();
};

//# sourceMappingURL=view-glycan-composition-peak-group-database-search.js.map
