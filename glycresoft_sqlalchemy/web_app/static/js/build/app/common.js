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
          $(".lean-overlay").remove();
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
