var ActionBook, DataSource, PartialSource, makeAPIGet, makeParameterizedAPIGet, makePartialGet;

ActionBook = {
  home: {
    container: '#home-layer',
    name: 'home-layer',
    closeable: false
  },
  addSample: {
    contentURL: '/add_sample',
    name: 'add-sample'
  },
  peakGroupingMatchSamples: {
    contentURL: '/peak_grouping_match_samples',
    name: "peak-grouping-match-samples"
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
    name: "view-database-search-results",
    method: "post"
  },
  viewHypothesis: {
    contentURLTemplate: "/view_hypothesis/{hypothesis_id}",
    method: "post"
  }
};

makeAPIGet = function(url) {
  return function(callback) {
    return $.get(url).success(callback);
  };
};

makeParameterizedAPIGet = function(url) {
  return function(params, callback) {
    return $.get(url.format(params)).success(callback);
  };
};

DataSource = {
  hypotheses: makeAPIGet("/api/hypotheses"),
  samples: makeAPIGet("/api/samples"),
  hypothesisSampleMatches: makeAPIGet("/api/hypothesis_sample_matches"),
  tasks: makeAPIGet("/api/tasks"),
  glycopeptideMatches: makeAPIGet("/api/glycopeptide_matches")
};

makePartialGet = function(url, method) {
  return function(parameters, callback) {
    return $[method](url.format(parameters)).success(callback);
  };
};

PartialSource = {
  glycopeptideCompositionDetailsModal: makePartialGet('/view_database_search_results/view_glycopeptide_composition_details/{id}', "get"),
  glycanCompositionDetailsModal: makePartialGet('/view_database_search_results/view_glycan_composition_details/{id}', "get")
};

//# sourceMappingURL=bind-urls.js.map

var Application, renderTask,
  extend = function(child, parent) { for (var key in parent) { if (hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; },
  hasProp = {}.hasOwnProperty;

Application = (function(superClass) {
  extend(Application, superClass);

  function Application(options) {
    var self;
    if (options == null) {
      options = {};
    }
    console.log("Instantiating Application", this);
    Application.__super__.constructor.call(this, options.actionContainer, options);
    this.version = [0, 0, 1];
    this.context = {};
    this.settings = {};
    this.tasks = {};
    this.sideNav = $('.side-nav');
    this.colors = new ColorManager();
    self = this;
    this.connectEventSource();
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
    this.handleMessage('task-error', (function(_this) {
      return function(data) {
        var task;
        task = self.tasks[data.id];
        task.status = 'error';
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
    this.on("layer-change", (function(_this) {
      return function(data) {
        return _this.colors.update();
      };
    })(this));
  }

  Application.prototype.connectEventSource = function() {
    console.log("Establishing EventSource connection");
    return this.eventStream = new EventSource('/stream');
  };

  Application.prototype.runInitializers = function() {
    var initializer, j, len, ref, results;
    ref = Application.initializers;
    results = [];
    for (j = 0, len = ref.length; j < len; j++) {
      initializer = ref[j];
      results.push(initializer.apply(this, null));
    }
    return results;
  };

  Application.prototype.updateSettings = function(payload) {
    if (payload == null) {
      payload = {};
    }
    return $.post('/preferences', payload).success((function(_this) {
      return function(data) {
        var k, v;
        for (k in data) {
          v = data[k];
          _this.settings[k] = v;
        }
        return _this.emit("update_settings");
      };
    })(this)).error(function(err) {
      return console.log("error in updateSettings", err, arguments);
    });
  };

  Application.prototype.updateTaskList = function() {
    var clickTask, doubleClickTask, self, taskListContainer;
    taskListContainer = this.sideNav.find('.task-list-container ul');
    clickTask = function(event) {
      var handle, id, state;
      handle = $(this);
      state = handle.attr('data-status');
      id = handle.attr('data-id');
      if (state === 'finished') {
        delete self.tasks[id];
        handle.fadeOut();
        handle.remove();
      }
    };
    self = this;
    doubleClickTask = function(event) {
      var handle, id;
      handle = $(this);
      id = handle.attr('data-id');
      return $.get("/internal/log/" + id, (function(_this) {
        return function(message) {
          return self.displayMessageModal(message);
        };
      })(this));
    };
    taskListContainer.html(_.map(this.tasks, renderTask).join(''));
    taskListContainer.find('li').map(function(i, li) {
      return contextMenu(li, {
        "View Log": doubleClickTask
      });
    });
    taskListContainer.find('li').click(clickTask);
    return taskListContainer.find("li").dblclick(doubleClickTask);
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
      var self;
      self = this;
      return $(function() {
        self.container = $(self.options.actionContainer);
        self.sideNav = $('.side-nav');
        self.addLayer(ActionBook.home);
        $("#run-matching").click(function(event) {
          return setupAjaxForm("/ms1_or_ms2_choice?ms1_choice=peakGroupingMatchSamples&ms2_choice=tandemMatchSamples", "#message-modal");
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
    }, function() {
      return this.handleMessage("files-to-download", (function(_this) {
        return function(data) {
          var file, j, len, ref, results;
          ref = data.files;
          results = [];
          for (j = 0, len = ref.length; j < len; j++) {
            file = ref[j];
            results.push(_this.downloadFile(file));
          }
          return results;
        };
      })(this));
    }, function() {
      return this.on("update_settings", (function(_this) {
        return function() {
          var layer;
          layer = _this.getShowingLayer();
          if (layer.name !== ActionBook.home.name) {
            console.log("Updated Settings, Current Layer:", layer.name);
            return layer.setup();
          }
        };
      })(this));
    }
  ];

  Application.prototype.loadData = function() {
    DataSource.hypotheses((function(_this) {
      return function(d) {
        console.log('hypothesis', d);
        _this.hypotheses = d;
        return _this.emit("render-hypotheses");
      };
    })(this));
    DataSource.samples((function(_this) {
      return function(d) {
        console.log('samples', d);
        _this.samples = d;
        return _this.emit("render-samples");
      };
    })(this));
    DataSource.hypothesisSampleMatches((function(_this) {
      return function(d) {
        console.log('hypothesisSampleMatches', d);
        _this.hypothesisSampleMatches = d;
        return _this.emit("render-hypothesis-sample-matches");
      };
    })(this));
    DataSource.tasks((function(_this) {
      return function(d) {
        console.log('tasks', d);
        _this.tasks = d;
        return _this.updateTaskList();
      };
    })(this));
    return this.colors.update();
  };

  Application.prototype.downloadFile = function(filePath) {
    return window.location = "/internal/file_download/" + btoa(filePath);
  };

  Application.prototype.displayMessageModal = function(message) {
    var container;
    container = $("#message-modal");
    container.find('.modal-content').html(message);
    return container.openModal();
  };

  Application.prototype.ajaxWithContext = function(url, options) {
    var data;
    if (options == null) {
      options = {
        data: {}
      };
    }
    data = options.data;
    data['settings'] = this.settings;
    data['context'] = this.context;
    options.method = "POST";
    options.data = JSON.stringify(data);
    options.contentType = "application/json";
    return $.ajax(url, options);
  };

  return Application;

})(ActionLayerManager);

renderTask = function(task) {
  return '<li data-id=\'{id}\' data-status=\'{status}\'><b>{name}</b> ({status})</li>'.format(task);
};

$(function() {
  var options;
  window.GlycReSoft = new Application(options = {
    actionContainer: ".action-layer-container"
  });
  console.log("updating Application");
  GlycReSoft.runInitializers();
  GlycReSoft.updateSettings();
  return GlycReSoft.updateTaskList();
});

//# sourceMappingURL=common.js.map

var ConstraintInputGrid, MonosaccharideInputWidgetGrid;

MonosaccharideInputWidgetGrid = (function() {
  MonosaccharideInputWidgetGrid.prototype.template = "<div class='monosaccharide-row row'>\n    <div class='input-field col s2'>\n        <label for='mass_shift_name'>Monosaccharide Name</label>\n        <input class='monosaccharide-name' type='text' name='monosaccharide_name' placeholder='Name'>\n    </div>\n    <div class='input-field col s2'>\n        <label for='monosaccharide_mass_delta'>Lower Bound</label>\n        <input class='lower-bound' type='number' name='monosaccharide_lower_bound' placeholder='Lower Bound'>\n    </div>\n    <div class='input-field col s2'>\n        <label for='monosaccharide_max_count'>Upper Bound</label>    \n        <input class='upper-bound' type='number' min='0' placeholder='Upper Bound' name='monosaccharide_upper_bound'>\n    </div>\n    <div class='input-field col s2'>\n        <label for='monosaccharide_composition'>Monosaccharide Composition</label>\n        <input class='monosaccharide-composition' type='text' name='monosaccharide_composition' placeholder='Composition'>\n    </div>\n</div>";

  function MonosaccharideInputWidgetGrid(container) {
    this.counter = 0;
    this.container = $(container);
    this.monosaccharides = {};
  }

  MonosaccharideInputWidgetGrid.prototype.update = function() {
    var entry, i, len, monosaccharides, notif, notify, ref, row;
    monosaccharides = {};
    ref = this.container.find(".monosaccharide-row");
    for (i = 0, len = ref.length; i < len; i++) {
      row = ref[i];
      row = $(row);
      console.log(row);
      entry = {
        name: row.find(".monosaccharide-name").val(),
        lower_bound: row.find(".lower-bound").val(),
        upper_bound: row.find(".upper-bound").val(),
        composition: row.find(".monosaccharide-composition").val()
      };
      if (entry.name === "") {
        continue;
      }
      if (entry.name in monosaccharides) {
        row.addClass("warning");
        notify = new TinyNotification(0, 0, "This monosaccharide is already present.", row);
        row.data("tinyNotification", notify);
        console.log(notify);
      } else {
        row.removeClass("warning");
        if (row.data("tinyNotification") != null) {
          notif = row.data("tinyNotification");
          notif.dismiss();
          row.data("tinyNotification", void 0);
        }
        monosaccharides[entry.name] = entry;
      }
    }
    console.log(monosaccharides);
    return this.monosaccharides = monosaccharides;
  };

  MonosaccharideInputWidgetGrid.prototype.addEmptyRowOnEdit = function(addHeader) {
    var callback, row, self;
    if (addHeader == null) {
      addHeader = false;
    }
    row = $(this.template);
    if (!addHeader) {
      row.find("label").remove();
    }
    this.container.append(row);
    row.data("counter", ++this.counter);
    self = this;
    callback = function(event) {
      if (row.data("counter") === self.counter) {
        self.addEmptyRowOnEdit(false);
      }
      return $(this).parent().find("label").removeClass("active");
    };
    row.find("input").change(callback);
    return row.find("input").change((function(_this) {
      return function() {
        return _this.update();
      };
    })(this));
  };

  MonosaccharideInputWidgetGrid.prototype.addRow = function(name, lower, upper, composition, addHeader) {
    var row;
    if (addHeader == null) {
      addHeader = false;
    }
    row = $(this.template);
    if (!addHeader) {
      row.find("label").remove();
    }
    this.counter += 1;
    row.find(".monosaccharide-name").val(name);
    row.find(".lower-bound").val(lower);
    row.find(".upper-bound").val(upper);
    row.find(".monosaccharide-composition").val(composition);
    this.container.append(row);
    row.find("input").change((function(_this) {
      return function() {
        return _this.update();
      };
    })(this));
    console.log(row);
    return this.update();
  };

  return MonosaccharideInputWidgetGrid;

})();

ConstraintInputGrid = (function() {
  ConstraintInputGrid.prototype.template = "<div class=\"monosaccharide-constraints-row row\">\n    <div class='input-field col s2'>\n        <label for='left_hand_side'>Name</label>\n        <input class='monosaccharide-name' type='text' name='left_hand_side' placeholder='Name'>\n    </div>\n    <div class='input-field col s1' style='padding-left: 2px;padding-right: 2px;'>\n        <select class='browser-default' name='operator'>\n            <option>=</option>\n            <option>!=</option>\n            <option>&gt;</option>\n            <option>&lt;</option>\n            <option>&gt;=</option>\n            <option>&lt;=</option>\n        </select>\n    </div>\n    <div class='input-field col s2'>\n        <label for='right_hand_side'>Name/Value</label>\n        <input class='monosaccharide-name' type='text' name='right_hand_side' placeholder='Name/Value'>\n    </div>\n</div>";

  function ConstraintInputGrid(container, monosaccharideGrid) {
    this.counter = 0;
    this.container = $(container);
    this.constraints = {};
    this.monosaccharideGrid = monosaccharideGrid;
  }

  ConstraintInputGrid.prototype.addEmptyRowOnEdit = function(addHeader) {
    var callback, row, self;
    if (addHeader == null) {
      addHeader = false;
    }
    row = $(this.template);
    if (!addHeader) {
      row.find("label").remove();
    }
    this.container.append(row);
    row.data("counter", ++this.counter);
    self = this;
    callback = function(event) {
      if (row.data("counter") === self.counter) {
        self.addEmptyRowOnEdit(false);
      }
      return $(this).parent().find("label").removeClass("active");
    };
    row.find("input").change(callback);
    return row.find("input").change((function(_this) {
      return function() {
        return _this.update();
      };
    })(this));
  };

  ConstraintInputGrid.prototype.update = function() {
    var constraints, entry, getMonosaccharide, i, len, notif, notify, ref, row;
    constraints = [];
    ref = this.container.find(".monosaccharide-constraints-row");
    for (i = 0, len = ref.length; i < len; i++) {
      row = ref[i];
      row = $(row);
      console.log(row);
      entry = {
        lhs: row.find("input[name='left_hand_side']").val(),
        operator: row.find("input[name='operator']").val(),
        rhs: row.find("input[name='right_hand_side']").val()
      };
      if (entry.lhs === "" || entry.rhs === "") {
        continue;
      }
      getMonosaccharide = function(name) {
        return /^(\d+)(.+)/.exec(name)[2];
      };
      if (!(getMonosaccharide(entry.lhs) in this.monosaccharideGrid.monosaccharides)) {
        row.addClass("warning");
        notify = new TinyNotification(0, 0, entry.lhs + " is not defined.", row);
        row.data("tinyNotification", notify);
        console.log(notify);
      } else if (!(getMonosaccharide(entry.rhs) in this.monosaccharideGrid.monosaccharides)) {
        row.addClass("warning");
        if (row.data("tinyNotification") != null) {
          notif = row.data("tinyNotification");
          notif.dismiss();
          row.data("tinyNotification", void 0);
        }
        notify = new TinyNotification(0, 0, entry.rhs + " is not defined.", row);
        row.data("tinyNotification", notify);
        console.log(notify);
      } else {
        row.removeClass("warning");
        if (row.data("tinyNotification") != null) {
          notif = row.data("tinyNotification");
          notif.dismiss();
          row.data("tinyNotification", void 0);
        }
      }
      constraints.push(entry);
    }
    console.log(constraints);
    return this.constraints = constraints;
  };

  return ConstraintInputGrid;

})();

//# sourceMappingURL=glycan-composition-builder-ui.js.map

Application.prototype.renderHypothesisSampleMatchListAt = function(container) {
  var chunks, hsm, row, self, template;
  chunks = [];
  template = (function() {
    var i, len, ref, results;
    ref = _.sortBy(_.values(this.hypothesisSampleMatches), function(o) {
      return o.id;
    });
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      hsm = ref[i];
      hsm.name = hsm.name != null ? hsm.name : "HypothesisSampleMatch:" + hsm.target_hypothesis.name + "@" + hsm.sample_run_name;
      row = $("<div data-id=" + hsm.id + " class='list-item clearfix'> <span class='handle'>" + hsm.id + ". " + (hsm.name.replace('_', ' ')) + "</span> <small class='right' style='display:inherit'> " + (hsm.hypothesis_sample_match_type.replace('HypothesisSampleMatch', '')) + " <a class='remove-hsm mdi-content-clear'></a> </small> </div>");
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
        handle = $(this);
        return console.log("Removal of HypothesisSampleMatch is not implemented.");
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

Application.prototype.renderHypothesisListAt = function(container) {
  var chunks, hypothesis, i, len, ref, row, self, template;
  chunks = [];
  template = '';
  self = this;
  ref = _.sortBy(_.values(this.hypotheses), function(o) {
    return o.id;
  });
  for (i = 0, len = ref.length; i < len; i++) {
    hypothesis = ref[i];
    if (hypothesis.is_decoy) {
      continue;
    }
    row = $("<div data-id=" + hypothesis.id + " class='list-item clearfix'> <span class='handle'>" + hypothesis.id + ". " + (hypothesis.name.replace('_', ' ')) + "</span> <small class='right' style='display:inherit'> " + (hypothesis.hypothesis_type != null ? hypothesis.hypothesis_type : '-') + " <a class='remove-hypothesis mdi-content-clear'></a> </small> </div>");
    chunks.push(row);
    row.click(function(event) {
      var handle, hypothesisId, layer;
      handle = $(this);
      hypothesisId = handle.attr("data-id");
      self.addLayer(ActionBook.viewHypothesis, {
        "hypothesis_id": hypothesisId
      });
      layer = self.lastAdded;
      return self.setShowingLayer(layer);
    });
    row.find(".remove-hypothesis").click(function(event) {
      var handle;
      return handle = $(this);
    });
  }
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

var MassShiftInputWidget;

MassShiftInputWidget = (function() {
  var addEmptyRowOnEdit, counter, template;
  template = "<div class='mass-shift-row row'>\n    <div class='input-field col s3'>\n        <label for='mass_shift_name'>Mass Shift Name</label>\n        <input class='mass-shift-name' type='text' name='mass_shift_name' placeholder='Name'>\n    </div>\n    <div class='input-field col s3'>\n        <label for='mass_shift_mass_delta'>Mass &Delta;</label>\n        <input class='mass-delta' type='number' name='mass_shift_mass_delta' step=\"0.0001\" placeholder='Mass Shift'>\n    </div>\n    <div class='input-field col s3'>\n        <label for='mass_shift_max_count'>Maximum Count</label>    \n        <input class='max-count' type='number' min='0' placeholder='Maximum Count' name='mass_shift_max_count'>\n    </div>\n</div>";
  counter = 0;
  addEmptyRowOnEdit = function(container, addHeader) {
    var callback, row;
    if (addHeader == null) {
      addHeader = true;
    }
    container = $(container);
    if (addHeader) {
      row = $(template);
    } else {
      row = $(template);
      row.find("label").remove();
    }
    container.append(row);
    row.data("counter", ++counter);
    callback = function(event) {
      if (row.data("counter") === counter) {
        addEmptyRowOnEdit(container, false);
      }
      return $(this).parent().find("label").removeClass("active");
    };
    return row.find("input").change(callback);
  };
  return addEmptyRowOnEdit;
})();

//# sourceMappingURL=mass-shift-ui.js.map

var MonosaccharideFilter;

MonosaccharideFilter = (function() {
  function MonosaccharideFilter(parent, residueNames, rules) {
    if (rules == null) {
      if (GlycReSoft.settings.monosaccharide_filters == null) {
        GlycReSoft.settings.monosaccharide_filters = {};
      }
      rules = GlycReSoft.settings.monosaccharide_filters;
    }
    this.container = $("<div></div>").addClass("row");
    $(parent).append(this.container);
    this.residueNames = residueNames;
    this.rules = rules;
  }

  MonosaccharideFilter.prototype.makeFilterWidget = function(residue) {
    var rendered, rule, sanitizeName, self, template;
    rule = this.rules[residue];
    if (rule == null) {
      rule = {
        minimum: 0,
        maximum: 10,
        include: true
      };
      this.rules[residue] = rule;
    }
    residue.name = residue;
    residue.sanitizeName = sanitizeName = residue.replace(/[\(\),]/g, "_");
    template = "<span class=\"col s2 monosaccharide-filter\" data-name='" + residue + "'>\n    <p style='margin: 0px; margin-bottom: -10px;'>\n        <input type=\"checkbox\" id=\"" + sanitizeName + "_include\" name=\"" + sanitizeName + "_include\"/>\n        <label for=\"" + sanitizeName + "_include\"><b>" + residue + "</b></label>\n    </p>\n    <p style='margin-top: 0px; margin-bottom: 0px;'>\n        <input id=\"" + sanitizeName + "_min\" type=\"number\" placeholder=\"Minimum " + residue + "\" style='width: 45px;' min=\"0\"\n               value=\"" + rule.minimum + "\" max=\"" + rule.maximum + "\" name=\"" + sanitizeName + "_min\"/> : \n        <input id=\"" + sanitizeName + "_max\" type=\"number\" placeholder=\"Maximum " + residue + "\" style='width: 45px;' min=\"0\"\n               value=\"" + rule.maximum + "\" name=\"" + sanitizeName + "_max\"/>\n    </p>\n</span>";
    self = this;
    rendered = $(template);
    rendered.find("#" + sanitizeName + "_min").change(function() {
      rule.minimum = parseInt($(this).val());
      return self.changed();
    });
    rendered.find("#" + sanitizeName + "_max").change(function() {
      rule.maximum = parseInt($(this).val());
      return self.changed();
    });
    rendered.find("#" + sanitizeName + "_include").prop("checked", rule.include).click(function() {
      rule.include = $(this).prop("checked");
      return self.changed();
    });
    return rendered;
  };

  MonosaccharideFilter.prototype.render = function() {
    var i, len, ref, residue, results, widget;
    ref = this.residueNames;
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      residue = ref[i];
      widget = this.makeFilterWidget(residue);
      results.push(this.container.append(widget));
    }
    return results;
  };

  MonosaccharideFilter.prototype.changed = _.debounce((function() {
    return GlycReSoft.emit("update_settings");
  }), 1000);

  return MonosaccharideFilter;

})();

//# sourceMappingURL=monosaccharide-composition-filter.js.map

Application.prototype.renderSampleListAt = function(container) {
  var chunks, row, sample, template;
  chunks = [];
  template = (function() {
    var i, len, ref, results;
    ref = _.sortBy(_.values(this.samples), function(o) {
      return o.id;
    });
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      sample = ref[i];
      row = $("<div data-name=" + sample.name + " class='list-item clearfix'> <span class='handle'>" + (sample.name.replace('_', ' ')) + "</span> <small class='right' style='display:inherit'> " + sample.sample_type + " <a class='remove-sample mdi-content-clear'></a> </small> </div>");
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

var viewGlycopeptideHypothesis;

viewGlycopeptideHypothesis = function(hypothesisId) {
  var currentPage, displayTable, proteinContainer, proteinId, setup, setupGlycopeptideTablePageHandler, updateCompositionTablePage, updateProteinChoice;
  displayTable = void 0;
  currentPage = 1;
  proteinContainer = void 0;
  proteinId = void 0;
  setup = function() {
    proteinContainer = $("#protein-container");
    $('.protein-list-table tbody tr').click(updateProteinChoice);
    return updateProteinChoice.apply($('.protein-list-table tbody tr'));
  };
  setupGlycopeptideTablePageHandler = function(page) {
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
    url = "/view_glycopeptide_hypothesis/protein_view/" + proteinId;
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
    url = "/view_glycopeptide_hypothesis/protein_view/" + proteinId + "/" + page;
    console.log(url);
    return GlycReSoft.ajaxWithContext(url).success(function(doc) {
      currentPage = page;
      displayTable.html(doc);
      return setupGlycopeptideTablePageHandler(page);
    });
  };
  return setup();
};

//# sourceMappingURL=view-glycopeptide-hypothesis.js.map

var viewPeakGroupingDatabaseSearchResults;

viewPeakGroupingDatabaseSearchResults = function() {
  var currentPage, currentProtein, downloadCSV, glycopeptideDetailsModal, glycopeptideTable, setup, setupGlycopeptideCompositionTablePageHandlers, showGlycopeptideCompositionDetailsModal, unload, updateGlycopeptideCompositionTablePage, updateProteinChoice;
  glycopeptideDetailsModal = void 0;
  glycopeptideTable = void 0;
  currentPage = 1;
  currentProtein = void 0;
  setup = function() {
    $('.protein-match-table tbody tr').click(updateProteinChoice);
    updateProteinChoice.apply($('.protein-match-table tbody tr'));
    return $("#save-csv-file").click(downloadCSV);
  };
  setupGlycopeptideCompositionTablePageHandlers = function(page) {
    if (page == null) {
      page = 1;
    }
    $('.glycopeptide-match-row').click(function() {
      var textSelection;
      textSelection = window.getSelection();
      if (!textSelection.toString()) {
        return showGlycopeptideCompositionDetailsModal.apply(this);
      }
    });
    $(':not(.disabled) .next-page').click(function() {
      return updateGlycopeptideCompositionTablePage(page + 1);
    });
    $(':not(.disabled) .previous-page').click(function() {
      return updateGlycopeptideCompositionTablePage(page - 1);
    });
    return $('.pagination li :not(.active)').click(function() {
      var nextPage;
      nextPage = $(this).attr("data-index");
      if (nextPage != null) {
        nextPage = parseInt(nextPage);
        return updateGlycopeptideCompositionTablePage(nextPage);
      }
    });
  };
  updateGlycopeptideCompositionTablePage = function(page) {
    var url;
    if (page == null) {
      page = 1;
    }
    url = "/view_database_search_results/glycopeptide_matches_composition_table/" + currentProtein + "/" + page;
    console.log(url);
    return GlycReSoft.ajaxWithContext(url).success(function(doc) {
      currentPage = page;
      glycopeptideTable.html(doc);
      return setupGlycopeptideCompositionTablePageHandlers(page);
    });
  };
  updateProteinChoice = function() {
    var handle, id;
    handle = $(this);
    currentProtein = id = handle.attr('data-target');
    $("#chosen-protein-container").html("<div class=\"progress\"><div class=\"indeterminate\"></div></div>").fadeIn();
    return GlycReSoft.ajaxWithContext('/view_database_search_results/protein_composition_view/' + id).success(function(doc) {
      var tabs;
      $('#chosen-protein-container').hide();
      $('#chosen-protein-container').html(doc).fadeIn();
      tabs = $('ul.tabs');
      tabs.tabs();
      GlycReSoft.context["current_protein"] = id;
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
      glycopeptideDetailsModal = $('#peptide-detail-modal');
      glycopeptideTable = $("#glycopeptide-table");
      return updateGlycopeptideCompositionTablePage(1);
    }).error(function(error) {
      return console.log(arguments);
    });
  };
  showGlycopeptideCompositionDetailsModal = function() {
    var handle, id;
    handle = $(this);
    id = handle.attr('data-target');
    console.log("ID " + id + ", Open Modal");
    return PartialSource.glycopeptideCompositionDetailsModal({
      "id": id
    }, function(doc) {
      glycopeptideDetailsModal.find('.modal-content').html(doc);
      $(".lean-overlay").remove();
      console.log("Opening Modal");
      return glycopeptideDetailsModal.openModal();
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

//# sourceMappingURL=view-peak-grouping-database-search.js.map

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
