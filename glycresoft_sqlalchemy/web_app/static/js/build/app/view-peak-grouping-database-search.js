var viewPeakGroupingDatabaseSearchResults;

viewPeakGroupingDatabaseSearchResults = function() {
  var peptideDetailsModal, setup, showGlycopeptideCompositionDetailsModal, unload, updateProteinChoice;
  peptideDetailsModal = void 0;
  setup = function() {
    $('.protein-match-table tbody tr').click(updateProteinChoice);
    return updateProteinChoice.apply($('.protein-match-table tbody tr'));
  };
  updateProteinChoice = function() {
    var handle, id;
    handle = $(this);
    id = handle.attr('data-target');
    $("#chosen-protein-container").html("<div class=\"progress\"><div class=\"indeterminate\"></div></div>").fadeIn();
    return $.post('/view_database_search_results/protein_composition_view/' + id, GlycReSoft.context).success(function(doc) {
      var tabs;
      $('#chosen-protein-container').hide();
      $('#chosen-protein-container').html(doc).fadeIn();
      tabs = $('ul.tabs');
      tabs.tabs();
      if (GlycReSoft.context['protein-view-active-tab'] !== void 0) {
        console.log(GlycReSoft.context['protein-view-active-tab']);
        $('ul.tabs').tabs('select_tab', GlycReSoft.context['protein-view-active-tab']);
      } else {
        $('ul.tabs').tabs('select_tab', 'protein-overview');
      }
      $('.indicator').addClass('indigo');
      $('ul.tabs .tab a').click(function() {
        return GlycReSoft.context['protein-view-active-tab'] = $(this).attr('href').slice(1);
      });
      peptideDetailsModal = $('#peptide-detail-modal');
      return $('.glycopeptide-match-row').click(showGlycopeptideCompositionDetailsModal);
    }).error(function(error) {
      return console.log(arguments);
    });
  };
  showGlycopeptideCompositionDetailsModal = function() {
    var handle, id;
    handle = $(this);
    id = handle.attr('data-target');
    console.log(id);
    return PartialSource.glycopeptideCompositionDetailsModal({
      "id": id
    }, function(doc) {
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

//# sourceMappingURL=view-peak-grouping-database-search.js.map
