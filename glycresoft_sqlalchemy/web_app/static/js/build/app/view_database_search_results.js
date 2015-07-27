var ViewDatabaseSearchResults, getGlycopeptideMatchDetails, peptideDetailsModal, showGlycopeptideDetailsModal, updateProteinChoice,
  extend = function(child, parent) { for (var key in parent) { if (hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; },
  hasProp = {}.hasOwnProperty;

ViewDatabaseSearchResults = (function(superClass) {
  extend(ViewDatabaseSearchResults, superClass);

  function ViewDatabaseSearchResults() {
    $('.protein-match-table tbody tr').click(updateProteinChoice);
    updateProteinChoice.apply($('.protein-match-table tbody tr'), null);
    return;
  }

  ViewDatabaseSearchResults.prototype.initGlycopeptideTooltip = function() {
    $('svg .glycopeptide').customTooltip(this.glycopeptideTooltipCallback, 'protein-view-tooltip');
    $('svg .modification path').customTooltip(this.modificationTooltipCallback, 'protein-view-tooltip');
  };

  ViewDatabaseSearchResults.prototype.glycopeptideTooltipCallback = function(handle) {
    var template;
    template = '<div> <span><b>MS2 Score:</b> {ms2-score}</span><br> <span><b>q-value:</b> {q-value}</span><br> <b>{sequence}</b> </div>';
    return template.format({
      'sequence': handle.attr('data-sequence'),
      'ms2-score': handle.attr('data-ms2-score'),
      'q-value': handle.attr('data-q-value')
    });
  };

  ViewDatabaseSearchResults.prototype.modificationTooltipCallback = function(handle) {
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

  return ViewDatabaseSearchResults;

})(EventEmitter);

peptideDetailsModal = void 0;

updateProteinChoice = function() {
  var handle, id;
  handle = $(this);
  id = handle.attr('data-target');
  console.log(id);
  $('#chosen-protein-container').fadeOut();
  return $.get('/view_database_search_results/protein_view/' + id).success(function(doc) {
    var tabs;
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

//# sourceMappingURL=view_database_search_results.js.map
