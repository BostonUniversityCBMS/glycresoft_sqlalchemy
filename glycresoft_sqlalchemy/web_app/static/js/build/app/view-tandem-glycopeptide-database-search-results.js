var TandemGlycopeptideDatabaseSearchResultsController;

TandemGlycopeptideDatabaseSearchResultsController = (function() {
  var downloadCSV, getGlycopeptideMatchDetails, glycopeptideTooltipCallback, modificationTooltipCallback;
  glycopeptideTooltipCallback = function(handle) {
    var template;
    template = '<div><table>\n<tr><td style=\'padding:3px;\'><b>MS2 Score:</b> {ms2-score}</td><td style=\'padding:3px;\'><b>Mass:</b> {calculated-mass}</td></tr>\n<tr><td style=\'padding:3px;\'><b>q-value:</b> {q-value}</td><td style=\'padding:3px;\'><b>Spectrum Matches:</b> {spectra-count}</td></tr>\n</table>\n<span>{sequence}</span>\n</div>';
    return template.format({
      'sequence': new PeptideSequence(handle.attr('data-sequence')).format(GlycReSoft.colors),
      'ms2-score': parseFloat(handle.attr('data-ms2-score')).toFixed(4),
      'q-value': handle.attr('data-q-value'),
      "calculated-mass": parseFloat(handle.attr("data-calculated-mass")).toFixed(4),
      "spectra-count": handle.attr("data-spectra-count")
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
  getGlycopeptideMatchDetails = function(id, callback) {
    return $.get('/api/glycopeptide_match/' + id, callback);
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
  TandemGlycopeptideDatabaseSearchResultsController = (function() {
    function TandemGlycopeptideDatabaseSearchResultsController() {
      this.currentPage = 1;
      this.peptideDetailsModal = void 0;
      this.glycopeptideTable = void 0;
      this.currentProtein = void 0;
      this.setup();
    }

    TandemGlycopeptideDatabaseSearchResultsController.prototype.setup = function() {
      var handle, last_id, last_selector, updateProteinChoice;
      updateProteinChoice = this.updateProteinChoiceCallback();
      $('.protein-match-table tbody tr').click(updateProteinChoice);
      last_id = GlycReSoft.context['protein_id'];
      last_selector = '.protein-match-table tbody tr[data-target="' + last_id + '"]';
      handle = $(last_selector);
      if (handle.length !== 0) {
        updateProteinChoice.apply(handle);
      } else {
        updateProteinChoice.apply($($('.protein-match-table tbody tr')[0]));
      }
      $(".tooltipped").tooltip();
      return $("#save-csv-file").click(downloadCSV);
    };

    TandemGlycopeptideDatabaseSearchResultsController.prototype.setupGlycopeptideTablePageHandlers = function(page) {
      var self;
      if (page == null) {
        page = 1;
      }
      self = this;
      $('.glycopeptide-match-row').click(function() {
        var textSelection;
        textSelection = window.getSelection();
        if (!textSelection.toString()) {
          return self.showGlycopeptideDetailsModalCallback().apply(this);
        }
      });
      $(':not(.disabled) .next-page').click(function() {
        return self.updateGlycopeptideTablePage(page + 1);
      });
      $(':not(.disabled) .previous-page').click(function() {
        return self.updateGlycopeptideTablePage(page - 1);
      });
      return $('.pagination li :not(.active)').click(function() {
        var nextPage;
        nextPage = $(this).attr("data-index");
        if (nextPage != null) {
          nextPage = parseInt(nextPage);
          return self.updateGlycopeptideTablePage(nextPage);
        }
      });
    };

    TandemGlycopeptideDatabaseSearchResultsController.prototype.updateGlycopeptideTablePage = function(page) {
      var url;
      if (page == null) {
        page = 1;
      }
      url = "/view_database_search_results/glycopeptide_match_table/" + this.currentProtein + "/" + page;
      return GlycReSoft.ajaxWithContext(url).success((function(_this) {
        return function(doc) {
          _this.currentPage = page;
          _this.glycopeptideTable.html(doc);
          return _this.setupGlycopeptideTablePageHandlers(page);
        };
      })(this));
    };

    TandemGlycopeptideDatabaseSearchResultsController.prototype.initGlycopeptideOverviewPlot = function() {
      var glycopeptide, self;
      glycopeptide = $('svg .glycopeptide');
      glycopeptide.customTooltip(glycopeptideTooltipCallback, 'protein-view-tooltip');
      self = this;
      glycopeptide.hover(function(event) {
        var baseColor, handle, newColor;
        handle = $(this);
        baseColor = handle.find("path").css("fill");
        newColor = '#74DEC5';
        handle.data("baseColor", baseColor);
        return handle.find("path").css("fill", newColor);
      }, function(event) {
        var handle;
        handle = $(this);
        return handle.find("path").css("fill", handle.data("baseColor"));
      });
      glycopeptide.click(function(event) {
        var handle, id;
        handle = $(this);
        id = handle.data("record-id");
        return $.get('/view_database_search_results/view_glycopeptide_details/' + id).success(function(doc) {
          self.peptideDetailsModal.find('.modal-content').html(doc);
          $(".lean-overlay").remove();
          return self.peptideDetailsModal.openModal();
        });
      });
      return $('svg .modification path').customTooltip(modificationTooltipCallback, 'protein-view-tooltip');
    };

    TandemGlycopeptideDatabaseSearchResultsController.prototype.updateProteinChoiceCallback = function() {
      var callback, self;
      self = this;
      return callback = function() {
        var handle, id;
        handle = $(this);
        $('.active-row').removeClass("active-row");
        handle.addClass("active-row");
        id = handle.attr('data-target');
        self.currentProtein = id;
        $('#chosen-protein-container').fadeOut();
        $("#loading-top-level-chosen-protein-container").fadeIn();
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
            $("#loading-top-level-chosen-protein-container").fadeOut();
            $('#chosen-protein-container').html(doc).fadeIn();
            tabs = $('ul.tabs');
            GlycReSoft.ajaxWithContext("/view_database_search_results/protein_view/" + id + "/protein_overview_panel").success(function(svg) {
              $("#protein-overview").html(svg);
              return self.initGlycopeptideOverviewPlot();
            });
            GlycReSoft.ajaxWithContext("/view_database_search_results/protein_view/" + id + "/microheterogeneity_plot_panel").success(function(svgGal) {
              return $("#site-distribution").html(svgGal);
            });
            tabs.tabs();
            if (GlycReSoft.context['protein-view-active-tab'] !== void 0) {
              $('ul.tabs').tabs('select_tab', GlycReSoft.context['protein-view-active-tab']);
            } else {
              $('ul.tabs').tabs('select_tab', 'protein-overview');
            }
            $('ul.tabs .tab a').click(function() {
              return GlycReSoft.context['protein-view-active-tab'] = $(this).attr('href').slice(1);
            });
            $('.indicator').addClass('indigo');
            self.glycopeptideTable = $("#glycopeptide-table");
            self.updateGlycopeptideTablePage(1);
            self.peptideDetailsModal = $('#peptide-detail-modal');
            return GlycReSoft.context['protein_id'] = id;
          },
          error: function(error) {
            return console.log(arguments);
          }
        });
      };
    };

    TandemGlycopeptideDatabaseSearchResultsController.prototype.showGlycopeptideDetailsModalCallback = function() {
      var callback, self;
      self = this;
      return callback = function() {
        var handle, id;
        handle = $(this);
        id = handle.attr('data-target');
        return $.get('/view_database_search_results/view_glycopeptide_details/' + id).success(function(doc) {
          self.peptideDetailsModal.find('.modal-content').html(doc);
          $(".lean-overlay").remove();
          return self.peptideDetailsModal.openModal();
        });
      };
    };

    return TandemGlycopeptideDatabaseSearchResultsController;

  })();
  return TandemGlycopeptideDatabaseSearchResultsController;
})();

//# sourceMappingURL=view-tandem-glycopeptide-database-search-results.js.map
