var ActionBook, DataSource, makeAPIGet;

ActionBook = {
  home: {
    container: '#home-layer',
    name: 'home-layer'
  },
  addSample: {
    contentURL: '/add_sample',
    name: 'add-sample'
  },
  tandemMatchSamples: {
    contentURL: '/tandem_match_samples',
    name: 'tandem-match-samples'
  },
  naiveGlycopeptideSearchSpace: {
    contentURL: "/glycopeptide_search_space",
    name: "glycopeptide-search-space"
  },
  naiveGlycanSearchSpace: {
    contentURL: "/glycan_search_space",
    name: "glycan-search-space"
  },
  viewDatabaseSearchResults: {
    contentURLTemplate: "/view_database_search_results/{hypothesis_sample_match_id}",
    name: "view-database-search-results"
  }
};

makeAPIGet = function(url) {
  return function(callback) {
    return $.get(url).success(callback);
  };
};

DataSource = {
  hypotheses: makeAPIGet("/api/hypotheses"),
  samples: makeAPIGet("/api/samples"),
  hypothesisSampleMatches: makeAPIGet("/api/hypothesis_sample_matches"),
  tasks: makeAPIGet("/api/tasks")
};

//# sourceMappingURL=bind-urls.js.map
;var Application, options, renderTask,
  extend = function(child, parent) { for (var key in parent) { if (hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; },
  hasProp = {}.hasOwnProperty;

Application = (function(superClass) {
  extend(Application, superClass);

  function Application(options) {
    var self;
    if (options == null) {
      options = {};
    }
    Application.__super__.constructor.call(this, options.actionContainer, options);
    this.version = [0, 0, 1];
    this.context = {};
    this.settings = {};
    this.tasks = {};
    this.sideNav = $('.side-nav');
    self = this;
    this.eventStream = new EventSource('/stream');
    this.handleMessage('update', (function(_this) {
      return function(data) {
        Materialize.toast(data.replace(/"/g, ''), 4000);
      };
    })(this));
    this.handleMessage('task-queued', (function(_this) {
      return function(data) {
        self.tasks[data.id] = {
          'id': data.id,
          'name': data.name,
          'status': 'queued'
        };
        self.updateTaskList();
      };
    })(this));
    this.handleMessage('task-start', (function(_this) {
      return function(data) {
        self.tasks[data.id] = {
          'id': data.id,
          'name': data.name,
          'status': 'running'
        };
        self.updateTaskList();
      };
    })(this));
    this.handleMessage('task-complete', (function(_this) {
      return function(data) {
        var err;
        try {
          self.tasks[data.id].status = 'finished';
        } catch (_error) {
          err = _error;
          self.tasks[data.id] = {
            'id': data.id,
            'name': data.name,
            'status': 'finished'
          };
        }
        self.updateTaskList();
      };
    })(this));
    this.handleMessage('new-sample', (function(_this) {
      return function(data) {
        _this.samples[data.id] = data;
        return _this.emit("render-samples");
      };
    })(this));
    this.handleMessage('new-hypothesis', (function(_this) {
      return function(data) {
        _this.hypotheses[data.id] = data;
        return _this.emit("render-hypotheses");
      };
    })(this));
    this.handleMessage('new-hypothesis-sample-match', (function(_this) {
      return function(data) {
        _this.hypothesisSampleMatches[data.id] = data;
        return _this.emit("render-hypothesis-sample-matches");
      };
    })(this));
  }

  Application.prototype.runInitializers = function() {
    var i, initializer, len, ref, results;
    ref = Application.initializers;
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      initializer = ref[i];
      results.push(initializer.apply(this, null));
    }
    return results;
  };

  Application.prototype.updateSettings = function() {
    return $.post('/internal/update_settings', this.settings).success(function(data) {
      return this.settings = data;
    }).error(function(err) {
      return console.log(err);
    });
  };

  Application.prototype.updateTaskList = function() {
    var clickTask, self, taskListContainer;
    taskListContainer = this.sideNav.find('.task-list-container ul');
    clickTask = function(event) {
      var handle, id, state;
      handle = $(this);
      state = handle.attr('data-status');
      id = handle.attr('data-id');
      if (state === 'finished') {
        console.log(self.tasks[id]);
        delete self.tasks[id];
        handle.fadeOut();
        handle.remove();
      }
    };
    taskListContainer.html(_.map(this.tasks, renderTask).join(''));
    self = this;
    return taskListContainer.find('li').click(clickTask);
  };

  Application.prototype.handleMessage = function(messageType, handler) {
    return this.eventStream.addEventListener(messageType, function(event) {
      var data;
      data = JSON.parse(event.data);
      return handler(data);
    });
  };

  Application.initializers = [
    function() {
      return console.log(this);
    }, function() {
      var self;
      self = this;
      return $(function() {
        self.container = $(self.options.actionContainer);
        self.sideNav = $('.side-nav');
        self.addLayer(ActionBook.home);
        $("#run-matching").click(function(event) {
          self.addLayer(ActionBook.tandemMatchSamples);
          return self.setShowingLayer(self.lastAdded);
        });
        $("#build-glycan-search-space").click(function(event) {
          self.addLayer(ActionBook.naiveGlycanSearchSpace);
          return self.setShowingLayer(self.lastAdded);
        });
        return $("#build-glycopeptide-search-space").click(function(event) {
          self.addLayer(ActionBook.naiveGlycopeptideSearchSpace);
          return self.setShowingLayer(self.lastAdded);
        });
      });
    }, function() {
      return this.loadData();
    }
  ];

  Application.prototype.loadData = function() {
    DataSource.hypotheses((function(_this) {
      return function(d) {
        _this.hypotheses = d;
        return _this.emit("render-hypotheses");
      };
    })(this));
    DataSource.samples((function(_this) {
      return function(d) {
        _this.samples = d;
        return _this.emit("render-samples");
      };
    })(this));
    return DataSource.hypothesisSampleMatches((function(_this) {
      return function(d) {
        _this.hypothesisSampleMatches = d;
        return _this.emit("render-hypothesis-sample-matches");
      };
    })(this));
  };

  return Application;

})(ActionLayerManager);

renderTask = function(task) {
  return '<li data-id=\'{id}\' data-status=\'{status}\'><b>{name}</b> ({status})</li>'.format(task);
};

window.GlycReSoft = new Application(options = {
  actionContainer: ".action-layer-container"
});

$(function() {
  console.log("updating Application");
  GlycReSoft.runInitializers();
  GlycReSoft.updateSettings();
  return GlycReSoft.updateTaskList();
});

//# sourceMappingURL=common.js.map
;Application.prototype.renderHypothesisSampleMatchListAt = function(container) {
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
      hsm.name = hsm.name != null ? hsm.name : "HypothesisSampleMatch:" + hsm.target_hypothesis.name + "@" + hsm.sample_run_name;
      row = $("<div data-id=" + hsm.id + " class=''> <span class='handle'>" + (hsm.name.replace('_', ' ')) + "</span> <small class='right'>" + " <a class='remove-hsm mdi-content-clear'></a></small></div>");
      chunks.push(row);
      self = this;
      row.click(function(event) {
        var handle, id;
        handle = $(this);
        id = handle.attr('data-id');
        self.addLayer(ActionBook.viewDatabaseSearchResults, {
          hypothesis_sample_match_id: id
        });
        console.log(self.layers);
        console.log(self.lastAdded);
        self.context["hypothesis_sample_match_id"] = id;
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
;Application.prototype.renderHypothesisListAt = function(container) {
  var chunks, hypothesis, row, template;
  chunks = [];
  template = (function() {
    var i, len, ref, results;
    ref = _.sortBy(_.values(this.hypotheses), function(o) {
      return o.name;
    });
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      hypothesis = ref[i];
      row = $("<div data-id=" + hypothesis.id + " class=''> <span class='handle'>" + (hypothesis.name.replace('_', ' ')) + "</span> <small class='right'>" + (hypothesis.hypothesis_type != null ? hypothesis.hypothesis_type : '-') + " <a class='remove-hypothesis mdi-content-clear'></a></small></div>");
      chunks.push(row);
      results.push(row.find(".remove-hypothesis").click(function(event) {
        var handle;
        return handle = $(this);
      }));
    }
    return results;
  }).call(this);
  return $(container).html(chunks);
};

Application.initializers.push(function() {
  return this.on("render-hypotheses", (function(_this) {
    return function() {
      return _this.renderHypothesisListAt(".hypothesis-list");
    };
  })(this));
});

//# sourceMappingURL=hypothesis-ui.js.map
;Application.prototype.renderSampleListAt = function(container) {
  var chunks, row, sample, template;
  chunks = [];
  template = (function() {
    var i, len, ref, results;
    ref = _.sortBy(_.values(this.samples), function(o) {
      return o.name;
    });
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      sample = ref[i];
      row = $("<div data-name=" + sample.name + "> <span class='handle'>" + (sample.name.replace('_', ' ')) + "</span> <small class='right'>" + sample.sample_type + " <a class='remove-sample mdi-content-clear'></a></small></div>");
      chunks.push(row);
      results.push(row.find(".remove-sample").click(function(event) {
        var handle;
        handle = $(this);
        return console.log(handle);
      }));
    }
    return results;
  }).call(this);
  return $(container).html(chunks);
};

Application.initializers.push(function() {
  return this.on("render-samples", (function(_this) {
    return function() {
      return _this.renderSampleListAt(".sample-list");
    };
  })(this));
});

//# sourceMappingURL=sample-ui.js.map
;var doZoom, viewDatabaseSearchResults;

doZoom = function() {
  var svg, zoom;
  svg = d3.select("svg g");
  zoom = function() {
    return svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
  };
  return d3.select("svg g").call(d3.behavior.zoom().scaleExtent([1, 8]).on("zoom", zoom));
};

viewDatabaseSearchResults = function() {
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
