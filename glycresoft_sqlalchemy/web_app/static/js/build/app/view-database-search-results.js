var doZoom, viewTandemGlycopeptideDatabaseSearchResults;

doZoom = function() {
  var svg, zoom;
  svg = d3.select("svg g");
  zoom = function() {
    return svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
  };
  return d3.select("svg g").call(d3.behavior.zoom().scaleExtent([1, 8]).on("zoom", zoom));
};

viewTandemGlycopeptideDatabaseSearchResults = function() {
  var getGlycopeptideMatchDetails, glycopeptideTooltipCallback, initGlycopeptideOverviewPlot, modificationTooltipCallback, peptideDetailsModal, setup, showGlycopeptideDetailsModal, unload, updateProteinChoice;
  peptideDetailsModal = void 0;
  setup = function() {
    $('.protein-match-table tbody tr').click(updateProteinChoice);
    return updateProteinChoice.apply($('.protein-match-table tbody tr'));
  };
  initGlycopeptideOverviewPlot = function() {
    $('svg .glycopeptide').customTooltip(glycopeptideTooltipCallback, 'protein-view-tooltip');
    return $('svg .modification path').customTooltip(modificationTooltipCallback, 'protein-view-tooltip');
  };
  glycopeptideTooltipCallback = function(handle) {
    var template;
    template = '<div> <span><b>MS2 Score:</b> {ms2-score}</span><br> <span><b>q-value:</b> {q-value}</span><br> <b>{sequence}</b> </div>';
    return template.format({
      'sequence': handle.attr('data-sequence'),
      'ms2-score': handle.attr('data-ms2-score'),
      'q-value': handle.attr('data-q-value')
    });
  };
  modificationTooltipCallback = function(handle) {
    var sequence, template, value;
    template = '<div> <span>{value}</span> </div>';
    value = handle.parent().attr('data-modification-type');
    if (value === 'HexNAc') {
      sequence = $('#' + handle.parent().attr('data-parent')).attr('data-sequence');
      value = 'HexNAc - Glycosylation: ' + sequence.split(/(\[|\{)/).slice(1).join('');
    }
    return template.format({
      'value': value
    });
  };
  updateProteinChoice = function() {
    var handle, id;
    handle = $(this);
    id = handle.attr('data-target');
    $("#chosen-protein-container").html("<div class=\"progress\"><div class=\"indeterminate\"></div></div>").fadeIn();
    return $.post('/view_database_search_results/protein_view/' + id, GlycReSoft.context).success(function(doc) {
      var tabs;
      $('#chosen-protein-container').hide();
      $('#chosen-protein-container').html(doc).fadeIn();
      initGlycopeptideOverviewPlot();
      tabs = $('ul.tabs');
      tabs.tabs();
      if (GlycReSoft.context['protein-view-active-tab'] !== void 0) {
        console.log(GlycReSoft.context['protein-view-active-tab']);
        $('ul.tabs').tabs('select_tab', GlycReSoft.context['protein-view-active-tab']);
      } else {
        $('ul.tabs').tabs('select_tab', 'protein-overview');
      }
      $('ul.tabs .tab a').click(function() {
        return GlycReSoft.context['protein-view-active-tab'] = $(this).attr('href').slice(1);
      });
      $('.indicator').addClass('indigo');
      $('.glycopeptide-match-row').click(showGlycopeptideDetailsModal);
      return peptideDetailsModal = $('#peptide-detail-modal');
    }).error(function(error) {
      return console.log(arguments);
    });
  };
  getGlycopeptideMatchDetails = function(id, callback) {
    return $.get('/api/glycopeptide_match/' + id, callback);
  };
  showGlycopeptideDetailsModal = function() {
    var handle, id;
    handle = $(this);
    id = handle.attr('data-target');
    return $.get('/view_database_search_results/view_glycopeptide_details/' + id).success(function(doc) {
      peptideDetailsModal.find('.modal-content').html(doc);
      $(".lean-overlay").remove();
      return peptideDetailsModal.openModal();
    });
  };
  unload = function() {
    return GlycReSoft.removeCurrentLayer();
  };
  return setup();
};

//# sourceMappingURL=view-database-search-results.js.map
