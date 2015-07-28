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
    this.handleMessage("new-sample", (function(_this) {
      return function(data) {
        return console.log(data);
      };
    })(this));
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
        console.log(self.options.actionContainer);
        self.sideNav = $('.side-nav');
        return self.addLayer(ActionBook.home);
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
;var SampleManager;

SampleManager = (function() {
  function SampleManager(samples, listingContainer) {
    this.samples = samples;
    this.container = listingContainer;
  }

  SampleManager.prototype.renderSampleList = function() {
    var chunks, i, len, ref, results, sample, template;
    chunks = [];
    template = "<tr> <td data-id={id}>{name} <small>{sample_type}</small></td> </tr>";
    ref = this.samples;
    results = [];
    for (i = 0, len = ref.length; i < len; i++) {
      sample = ref[i];
      results.push(chunks.push(template.format(sample)));
    }
    return results;
  };

  return SampleManager;

})();

//# sourceMappingURL=sample_management.js.map
;var viewDatabaseSearchResults;

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
