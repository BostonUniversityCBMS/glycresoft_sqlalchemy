var ActionBook, ActionLayer, ActionLayerManager,
  extend = function(child, parent) { for (var key in parent) { if (hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; },
  hasProp = {}.hasOwnProperty;

ActionLayerManager = (function(superClass) {
  extend(ActionLayerManager, superClass);

  ActionLayerManager.HOME_LAYER = "home-layer";

  function ActionLayerManager(container, options) {
    this.container = $(container);
    this.options = options;
    this.layers = {};
    this.layerCounter = 0;
    this.on('layer-change', function(event) {
      materialRefresh();
    });
    return;
  }

  ActionLayerManager.prototype.incLayerCounter = function() {
    this.layerCounter++;
    return this.layerCounter;
  };

  ActionLayerManager.prototype.add = function(layer) {
    this.layers[layer.id] = layer;
    this.container.append(layer.container);
  };

  ActionLayerManager.prototype.get = function(id) {
    return this.layers[id];
  };

  ActionLayerManager.prototype.getShowingLayer = function() {
    var result;
    result = void 0;
    _.forIn(this.layers, function(layer, id) {
      if (layer.showing) {
        result = layer;
      }
    });
    return result;
  };

  ActionLayerManager.prototype.setShowingLayer = function(id) {
    var current, next;
    current = this.getShowingLayer();
    next = this.get(id);
    try {
      current.hide();
    } catch (_error) {}
    try {
      next.show();
    } catch (_error) {
      this.get(ActionLayerManager.HOME_LAYER).show();
    }
    return this.emit('layer-change', {
      'layer': next
    });
  };

  ActionLayerManager.prototype.addLayer = function(options, params) {
    return new ActionLayer(this, options, params);
  };

  ActionLayerManager.prototype.removeLayer = function(id) {
    this.layers[id].container.remove();
    return delete this.layers[id];
  };

  ActionLayerManager.prototype.removeCurrentLayer = function(next) {
    var current;
    if (next == null) {
      next = ActionLayerManager.HOME_LAYER;
    }
    current = this.getShowingLayer();
    this.setShowingLayer(next);
    return current.dispose();
  };

  return ActionLayerManager;

})(EventEmitter);

ActionLayer = (function() {
  ActionLayer.actions = {};

  function ActionLayer(manager, options, params) {
    this.manager = manager;
    this.options = options;
    this.params = params;
    this.contentURL = options.contentURL;
    if (!options.container) {
      if (this.params != null) {
        this.id = options.name + "-" + manager.incLayerCounter();
      } else {
        this.id = options.name || 'action-layer-' + manager.incLayerCounter();
      }
      this.container = $('<div></div>').attr('id', this.id);
      this.setup();
    } else {
      this.container = $(options.container);
      this.id = this.container.attr('id');
    }
    this.name = options.name || 'layer-' + this.id;
    this.container.attr('data-name', this.name).addClass('container');
    this.manager.add(this);
    this.showing = false;
    this.hide();
    if (this.manager.getShowingLayer() === void 0) {
      this.show();
    }
  }

  ActionLayer.prototype.setup = function() {
    var self;
    console.log("Setting up", this);
    if (this.options.contentURLTemplate != null) {
      this.contentURL = this.options.contentURLTemplate.format(this.params);
    }
    self = this;
    return $.get(this.contentURL).success(function(doc) {
      self.container.hide();
      self.container.html(doc);
      return self.container.find('script').each(function(i, tag) {
        var srcURL;
        tag = $(tag);
        srcURL = tag.attr('src');
        console.log("Setting up script", tag);
        if (srcURL !== void 0) {
          return $.getScript(srcURL);
        }
      });
    });
  };

  ActionLayer.prototype.show = function() {
    this.container.fadeIn(100);
    this.showing = true;
  };

  ActionLayer.prototype.hide = function() {
    this.container.fadeOut(100);
    this.showing = false;
  };

  ActionLayer.prototype.dispose = function() {
    this.container.remove();
    return delete this.manager.layers[this.id];
  };

  return ActionLayer;

})();

ActionBook = {
  home: {
    container: '#home-layer',
    name: 'home-layer'
  },
  addSample: {
    contentURL: '/add_sample',
    name: 'add-sample'
  },
  matchSamples: {
    contentURL: '/match_samples',
    name: 'match-samples'
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
    contentURLTemplate: "/view_database_search_results/{}",
    name: "view-database-search-results"
  }
};

//# sourceMappingURL=action-layer.js.map
;var ajaxForm, setupAjaxForm;

ajaxForm = function(formHandle, success, error) {
  return $(formHandle).on('submit', function(event) {
    var ajaxParams, data, encoding, handle, method, url;
    event.preventDefault();
    handle = $(this);
    url = handle.attr('action');
    method = handle.attr('method');
    data = new FormData(this);
    encoding = handle.attr('enctype') || 'application/x-www-form-urlencoded; charset=UTF-8';
    ajaxParams = {
      'url': url,
      'method': method,
      'data': data,
      'processData': false,
      'contentType': false,
      'success': success,
      'error': error
    };
    return $.ajax(ajaxParams);
  });
};

setupAjaxForm = function(sourceUrl, container) {
  var isModal;
  container = $(container);
  isModal = container.hasClass('modal');
  $.get(sourceUrl).success(function(doc) {
    if (isModal) {
      container.find('.modal-content').html(doc);
      container.openModal();
      return container.find('form').submit(function(event) {
        return container.closeModal();
      });
    } else {
      return container.html(doc);
    }
  });
  return container.find('script').each(function(i, tag) {
    var srcURL;
    tag = $(tag);
    srcURL = tag.attr('src');
    if (srcURL !== void 0) {
      return $.getScript(srcURL);
    } else {
      return eval(tag.text());
    }
  });
};

//# sourceMappingURL=ajax-form.js.map
;$(function() {
  var $body, $tooltip, closeTooltip, openTooltip;
  $body = $('body');
  $tooltip = $('<div></div>').hide().css({
    'position': 'absolute',
    'z-index': '10'
  });
  openTooltip = function(event) {
    var content, handle;
    handle = $(this);
    content = handle.data('tooltip-content');
    if (typeof content === 'function') {
      content = content(handle);
    }
    content = content === void 0 ? 'This is a simple tooltip' : content;
    $tooltip.html(content).addClass(handle.data('tooltip-css-class')).css('top', event.pageY - 10 + 'px').css('left', event.pageX + 10 + 'px').show();
  };
  closeTooltip = function(event) {
    var handle;
    handle = $(this);
    $tooltip.html('').removeClass(handle.data('tooltip-css-class')).hide();
  };
  $body.append($tooltip);
  jQuery.fn.customTooltip = function(content, cssClass) {
    var handle;
    handle = $(this);
    if (content !== void 0) {
      handle.data('tooltip-content', content);
    }
    if (cssClass !== void 0) {
      handle.data('tooltip-css-class', cssClass);
    }
    handle.hover(openTooltip, closeTooltip);
  };
});

//# sourceMappingURL=custom-tooltip.js.map
;(function() {

  /*
  Implements {named} replacements similar to the simple format() method of strings from Python
   */
  String.prototype.format = function() {
    var data, i, keys, res;
    data = arguments;
    i = 0;
    keys = Object.keys(arguments);
    if (arguments.length === 1 && typeof arguments[0] === 'object') {
      data = arguments[0];
      keys = Object.keys(arguments);
    }
    res = this.replace(/\{([^\}]*)\}/g, function(placeholder, name, position) {
      var err, v;
      if (name === '') {
        name = keys[i];
        i++;
      }
      try {
        v = JSON.stringify(data[name]);
        if (v.length > 1) {
          v = v.slice(1, -1);
        }
        return v;
      } catch (_error) {
        err = _error;
        console.log(err, name, data);
        return void 0;
      }
    });
    return res;
  };
})();

//# sourceMappingURL=formatstring.js.map
;
/*
The file name mirroring done by the .file-field group in Materialize is set up on page load.
When these elements are added dynamically, they must be configured manually.

This code is taken from https://github.com/Dogfalo/materialize/blob/master/js/forms.js#L156
 */
var materialFileInput, materialRefresh;

materialRefresh = function() {
  try {
    materialFileInput();
  } catch (_error) {}
  try {
    Materialize.updateTextFields();
  } catch (_error) {}
};

materialFileInput = function() {
  $(document).on('change', '.file-field input[type="file"]', function() {
    var file_field, file_names, files, i, path_input;
    file_field = $(this).closest('.file-field');
    path_input = file_field.find('input.file-path');
    files = $(this)[0].files;
    file_names = [];
    i = 0;
    while (i < files.length) {
      file_names.push(files[i].name);
      i++;
    }
    path_input.val(file_names.join(', '));
    path_input.trigger('change');
  });
};

//# sourceMappingURL=material-shim.js.map
