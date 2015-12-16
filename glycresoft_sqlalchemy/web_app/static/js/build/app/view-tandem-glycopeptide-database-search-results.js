var viewTandemGlycopeptideDatabaseSearchResults;

viewTandemGlycopeptideDatabaseSearchResults = function() {
  var downloadCSV, getGlycopeptideMatchDetails, glycopeptideTooltipCallback, initGlycopeptideOverviewPlot, modificationTooltipCallback, peptideDetailsModal, setup, showGlycopeptideDetailsModal, updateProteinChoice;
  peptideDetailsModal = void 0;
  setup = function() {
    var handle, last_id, last_selector;
    $('.protein-match-table tbody tr').click(updateProteinChoice);
    last_id = GlycReSoft.context['protein_id'];
    last_selector = '.protein-match-table tbody tr[data-target="' + last_id + '"]';
    console.log(last_selector);
    handle = $(last_selector);
    console.log(handle);
    if (handle.length !== 0) {
      updateProteinChoice.apply(handle);
    } else {
      updateProteinChoice.apply($($('.protein-match-table tbody tr')[0]));
    }
    $(".tooltipped").tooltip();
    return $("#save-csv-file").click(downloadCSV);
  };
  initGlycopeptideOverviewPlot = function() {
    var glycopeptide;
    glycopeptide = $('svg .glycopeptide');
    glycopeptide.customTooltip(glycopeptideTooltipCallback, 'protein-view-tooltip');
    glycopeptide.click(function(event) {
      var handle, id;
      handle = $(this);
      id = handle.data("record-id");
      return $.get('/view_database_search_results/view_glycopeptide_details/' + id).success(function(doc) {
        peptideDetailsModal.find('.modal-content').html(doc);
        $(".lean-overlay").remove();
        return peptideDetailsModal.openModal();
      });
    });
    return $('svg .modification path').customTooltip(modificationTooltipCallback, 'protein-view-tooltip');
  };
  glycopeptideTooltipCallback = function(handle) {
    var template;
    template = '<div> <span><b>MS2 Score:</b> {ms2-score}</span><br> <span><b>q-value:</b> {q-value}</span><br> <span>{sequence}</span> </div>';
    return template.format({
      'sequence': new PeptideSequence(handle.attr('data-sequence')).format(GlycReSoft.colors),
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
    $('.active-row').removeClass("active-row");
    handle.addClass("active-row");
    id = handle.attr('data-target');
    $("#chosen-protein-container").html("<div class=\"progress\"><div class=\"indeterminate\"></div></div>").fadeIn();
    return $.ajax('/view_database_search_results/protein_view/' + id, {
      data: JSON.stringify({
        "context": GlycReSoft.context,
        "settings": GlycReSoft.settings
      }),
      contentType: "application/json",
      type: 'POST',
      success: function(doc) {
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
        $('.glycopeptide-match-row').click(function() {
          var textSelection;
          textSelection = window.getSelection();
          if (!textSelection.toString()) {
            return showGlycopeptideDetailsModal.apply(this);
          }
        });
        peptideDetailsModal = $('#peptide-detail-modal');
        return GlycReSoft.context['protein_id'] = id;
      },
      error: function(error) {
        return console.log(arguments);
      }
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

//# sourceMappingURL=view-tandem-glycopeptide-database-search-results.js.map
