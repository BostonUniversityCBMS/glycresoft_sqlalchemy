var viewDatabaseSearchResults;

viewDatabaseSearchResults = function() {
  var getGlycopeptideMatchDetails, glycopeptideTooltipCallback, initGlycopeptideTooltip, modificationTooltipCallback, peptideDetailsModal, setup, showGlycopeptideDetailsModal, unload, updateProteinChoice;
  peptideDetailsModal = void 0;
  setup = function() {
    $('.protein-match-table tbody tr').click(updateProteinChoice);
    return updateProteinChoice.apply($('.protein-match-table tbody tr'));
  };
  initGlycopeptideTooltip = function() {
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
    console.log(id, "About to fadeOut");
    $('#chosen-protein-container').fadeOut();
    return $.get('/view_database_search_results/protein_view/' + id).success(function(doc) {
      var tabs;
      console.log("About to fadeIn", $('#chosen-protein-container'));
      $('#chosen-protein-container').html(doc).fadeIn();
      initGlycopeptideTooltip();
      tabs = $('ul.tabs');
      tabs.tabs();
      if (GlycReSoft.context['protein-view-active-tab'] !== void 0) {
        console.log(GlycReSoft.context['protein-view-active-tab']);
        $('ul.tabs').tabs('select_tab', GlycReSoft.context['protein-view-active-tab']);
      } else {
        $('ul.tabs').tabs('select_tab', 'protein-overview');
      }
      $('ul.tabs .tab a').click(function() {
        GlycReSoft.context['protein-view-active-tab'] = $(this).attr('href').slice(1);
      });
      $('.indicator').addClass('indigo');
      $('.glycopeptide-match-row').click(showGlycopeptideDetailsModal);
      return peptideDetailsModal = $('#peptide-detail-modal');
    }).error(function(error) {
      return console.log(arguments);
    });
  };
  getGlycopeptideMatchDetails = function(id, callback) {
    $.get('/api/glycopeptide_match/' + id, callback);
  };
  showGlycopeptideDetailsModal = function() {
    var handle, id;
    handle = $(this);
    id = handle.attr('data-target');
    return $.get('/view_database_search_results/view_glycopeptide_details/' + id).success(function(doc) {
      peptideDetailsModal.find('.modal-content').html(doc);
      return peptideDetailsModal.openModal();
    });
  };
  unload = function() {
    return GlycReSoft.removeCurrentLayer();
  };
  return setup();
};

//# sourceMappingURL=view_database_search_results.js.map
