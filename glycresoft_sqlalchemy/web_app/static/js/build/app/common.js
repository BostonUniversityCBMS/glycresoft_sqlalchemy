var Application, options, renderTask,
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
