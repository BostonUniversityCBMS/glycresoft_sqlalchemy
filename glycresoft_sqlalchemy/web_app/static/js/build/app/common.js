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
    this.eventStream.addEventListener('update', function(event) {
      Materialize.toast(event.data.replace(/"/g, ''), 4000);
    });
    this.eventStream.addEventListener('task-queued', function(event) {
      var data;
      data = JSON.parse(event.data);
      self.tasks[data.id] = {
        'id': data.id,
        'name': data.name,
        'status': 'queued'
      };
      self.updateTaskList();
    });
    this.eventStream.addEventListener('task-start', function(event) {
      var data;
      data = JSON.parse(event.data);
      self.tasks[data.id] = {
        'id': data.id,
        'name': data.name,
        'status': 'running'
      };
      self.updateTaskList();
    });
    this.eventStream.addEventListener('task-complete', function(event) {
      var data, err;
      data = JSON.parse(event.data);
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
    });
  }

  Application.prototype.runInitializers = function() {
    var i, initializer, len, ref, results;
    console.log(Application, Application.initializers);
    ref = Application.initializers;
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      initializer = ref[i];
      console.log(initializer);
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

  Application.initializers = [
    function() {
      return console.log(this);
    }, function() {
      var self;
      self = this;
      return $(function() {
        self.container = $(self.options.actionContainer);
        console.log(self.options.actionContainer);
        self.sideNav = $('.side-nav');
        self.addLayer({
          'container': '#home-layer',
          'name': 'home-layer'
        });
        self.addLayer({
          'contentURL': '/add_sample',
          'name': 'add-sample'
        });
        return self.addLayer({
          'contentURL': '/match_samples',
          'name': 'match-samples'
        });
      });
    }
  ];

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
