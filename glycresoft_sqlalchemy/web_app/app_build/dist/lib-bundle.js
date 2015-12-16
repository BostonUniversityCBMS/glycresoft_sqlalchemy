var ActionLayer, ActionLayerManager,
  extend = function(child, parent) { for (var key in parent) { if (hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; },
  hasProp = {}.hasOwnProperty;

ActionLayerManager = (function(superClass) {
  extend(ActionLayerManager, superClass);

  ActionLayerManager.HOME_LAYER = "home-layer";

  function ActionLayerManager(container, options) {
    this.container = $(container);
    this.options = options;
    this.layers = {};
    this.lastAdded = null;
    this.layerCounter = 0;
    this.on('layer-change', function(event) {
      return materialRefresh();
    });
    this.layerStack = [];
  }

  ActionLayerManager.prototype.incLayerCounter = function() {
    this.layerCounter++;
    return this.layerCounter;
  };

  ActionLayerManager.prototype.add = function(layer) {
    if (layer.options.closeable == null) {
      layer.options.closeable = true;
    }
    this.layers[layer.id] = layer;
    this.container.append(layer.container);
    this.emit("layer-added", {
      "layer": layer
    });
    this.lastAdded = layer.id;
    return this;
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
    var current, i, next;
    console.log(id);
    current = this.getShowingLayer();
    next = this.get(id);
    try {
      current.hide();
    } catch (_error) {}
    try {
      next.show();
      i = this.findLayer(next);
      if (i !== -1) {
        this.layerStack.pop(i);
      }
      this.layerStack.push(next);
    } catch (_error) {
      this.get(ActionLayerManager.HOME_LAYER).show();
    }
    return this.emit('layer-change', {
      'layer': next
    });
  };

  ActionLayerManager.prototype.addLayer = function(options, params) {
    var layer;
    if (options.closeable == null) {
      options.closeable = true;
    }
    layer = new ActionLayer(this, options, params);
    console.log(layer);
    if (this.layerStack.length === 0) {
      this.layerStack.push(layer);
    }
    return this;
  };

  ActionLayerManager.prototype.removeLayer = function(id) {
    var i;
    this.layers[id].container.remove();
    i = this.findLayer(this.layers[id]);
    if (i !== -1) {
      this.layerStack.pop(i);
    }
    delete this.layers[id];
    return this;
  };

  ActionLayerManager.prototype.removeCurrentLayer = function(next) {
    var current;
    if (next == null) {
      next = null;
    }
    current = this.getShowingLayer();
    this.layerStack.pop();
    if (next == null) {
      next = this.layerStack[this.layerStack.length - 1];
    }
    this.setShowingLayer(next);
    current.dispose();
    return this;
  };

  ActionLayerManager.prototype.findLayer = function(targetLayer) {
    var index, j, layer, len, ref;
    index = -1;
    ref = this.layerStack;
    for (j = 0, len = ref.length; j < len; j++) {
      layer = ref[j];
      index += 1;
      if (layer.id === targetLayer.id) {
        return index;
      }
    }
    return -1;
  };

  return ActionLayerManager;

})(EventEmitter);

ActionLayer = (function() {
  ActionLayer.actions = {};

  function ActionLayer(manager, options, params, method) {
    if (method == null) {
      method = 'get';
    }
    this.manager = manager;
    this.options = options;
    this.params = params;
    this.contentURL = options.contentURL;
    this.method = method;
    if (this.options.method != null) {
      this.method = this.options.method;
    }
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
    var callback;
    console.log("Setting up", this);
    if (this.options.contentURLTemplate != null) {
      this.contentURL = this.options.contentURLTemplate.format(this.params);
    }
    callback = (function(_this) {
      return function(doc) {
        if (!_this.showing) {
          _this.container.hide();
        }
        _this.options.document = doc;
        _this.container.html(doc);
        _this.container.find('script').each(function(i, tag) {
          var srcURL;
          tag = $(tag);
          srcURL = tag.attr('src');
          console.log("Setting up script", tag);
          if (srcURL !== void 0) {
            return $.getScript(srcURL);
          }
        });
        materialRefresh();
        console.log("This layer can be closed? " + _this.options.closeable);
        if (_this.options.closeable) {
          return _this.container.prepend("<div>\n    <a class='dismiss-layer mdi-content-clear' onclick='GlycReSoft.removeCurrentLayer()'></a>\n</div>");
        }
      };
    })(this);
    if (this.method === "get") {
      return $.get(this.contentURL).success(callback);
    } else if (this.method === "post") {
      return $.ajax(this.contentURL, {
        contentType: "application/json",
        data: JSON.stringify({
          params: this.params,
          context: this.manager.context,
          settings: this.manager.settings
        }),
        success: callback,
        type: "POST"
      });
    }
  };

  ActionLayer.prototype.reload = function() {
    this.container.html(this.options.document);
    this.container.find('script').each(function(i, tag) {
      var srcURL;
      tag = $(tag);
      srcURL = tag.attr('src');
      console.log("Setting up script", tag);
      if (srcURL !== void 0) {
        return $.getScript(srcURL);
      }
    });
    materialRefresh();
    if (this.options.closeable) {
      return this.container.prepend("<div>\n<a class='dismiss-layer mdi-content-clear' onclick='GlycReSoft.removeCurrentLayer()'></a>\n</div>");
    }
  };

  ActionLayer.prototype.show = function() {
    this.container.fadeIn(100);
    return this.showing = true;
  };

  ActionLayer.prototype.hide = function() {
    this.container.fadeOut(100);
    return this.showing = false;
  };

  ActionLayer.prototype.dispose = function() {
    this.container.remove();
    return delete this.manager.layers[this.id];
  };

  return ActionLayer;

})();

//# sourceMappingURL=action-layer.js.map

var ajaxForm, setupAjaxForm;

ajaxForm = function(formHandle, success, error, transform) {
  console.log("Ajaxifying ", formHandle);
  return $(formHandle).on('submit', function(event) {
    var ajaxParams, data, encoding, handle, method, url;
    console.log(formHandle, "submitting...");
    event.preventDefault();
    handle = $(this);
    if (transform == null) {
      transform = function(form) {
        return new FormData(form);
      };
    }
    url = handle.attr('action');
    method = handle.attr('method');
    data = transform(this);
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

var contextMenu;

contextMenu = function(target, options, callback) {
  if (callback == null) {
    callback = null;
  }
  $(target).off("contextmenu", false);
  $(document).off("mousedown", false);
  return $(target).on("contextmenu", function(event) {
    var action, handle, item;
    event.preventDefault();
    handle = $(".context-menu");
    handle.empty();
    if (callback != null) {
      callback(handle);
    }
    for (item in options) {
      action = options[item];
      handle.append($("<li></li>").text(item).attr("data-action", item));
    }
    $(".context-menu li").click(function(e) {
      handle = $(this);
      action = options[handle.attr("data-action")];
      return action.apply(target);
    });
    return $(".context-menu").finish().toggle(100).css({
      top: event.pageY + 'px',
      left: event.pageX + 'px'
    });
  });
};

$(document).on("mousedown", function(e) {
  if (!$(e.target).parents(".context-menu").length > 0) {
    return $(".context-menu").hide(100);
  }
});

//# sourceMappingURL=context-menu.js.map

$(function() {
  var $body, $tooltip, closeTooltip, openTooltip, xOffset, yOffset;
  yOffset = -3;
  xOffset = 3;
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
    $tooltip.html(content).addClass(handle.data('tooltip-css-class')).css('top', event.pageY + yOffset + 'px').css('left', event.pageX + xOffset + 'px').show();
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

(function() {

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

var getProteinName, getProteinNamesFromMzIdentML, identifyProteomicsFormat;

identifyProteomicsFormat = function(file, callback) {
  var isMzidentML, reader;
  isMzidentML = function(lines) {
    var i, len, line;
    for (i = 0, len = lines.length; i < len; i++) {
      line = lines[i];
      if (/mzIdentML/.test(line)) {
        return true;
      }
    }
    return false;
  };
  reader = new FileReader();
  reader.onload = function() {
    var lines, proteomicsFileType;
    lines = this.result.split("\n");
    console.log(lines);
    proteomicsFileType = "fasta";
    if (isMzidentML(lines)) {
      proteomicsFileType = "mzIdentML";
    }
    return callback(proteomicsFileType);
  };
  return reader.readAsText(file.slice(0, 100));
};

getProteinName = function(line) {
  return line.split("_", 2)[1];
};

getProteinNamesFromMzIdentML = function(file, callback) {
  var chunksize, fr, offset, proteins, seek;
  fr = new FileReader();
  chunksize = 1024 * 8;
  offset = 0;
  proteins = {};
  fr.onload = function() {
    var i, len, line, lines, name;
    lines = this.result.split("\n");
    for (i = 0, len = lines.length; i < len; i++) {
      line = lines[i];
      if (/<ProteinDetectionHypothesis/i.test(line)) {
        name = getProteinName(line);
        console.log(name);
        proteins[name] = true;
      }
    }
    return seek();
  };
  fr.onerror = function(error) {
    return console.log(error);
  };
  seek = function() {
    if (offset >= file.size) {
      return callback(Object.keys(proteins));
    } else {
      fr.readAsText(file.slice(offset, offset + chunksize));
      return offset += chunksize / 2;
    }
  };
  return seek();
};

//# sourceMappingURL=infer-protein-data-format.js.map


/*
The file name mirroring done by the .file-field group in Materialize is set up on page load.
When these elements are added dynamically, they must be configured manually.

This code is taken from https://github.com/Dogfalo/materialize/blob/master/js/forms.js#L156
 */
var materialFileInput, materialRefresh;

materialRefresh = function() {
  try {
    $('select').material_select();
  } catch (_error) {}
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

var TinyNotification;

TinyNotification = (function() {
  TinyNotification.prototype.template = "<div class='notification-container'>\n    <a class='dismiss-notification mdi-content-clear'></a>\n    <div class='notification-content'>\n    </div>\n</div>";

  function TinyNotification(top, left, message, parent, css) {
    if (parent == null) {
      parent = 'body';
    }
    if (css == null) {
      css = {};
    }
    this.top = top;
    this.left = left;
    this.parent = parent;
    this.message = message;
    this.container = $(this.template);
    this.container.find(".notification-content").html(this.message);
    this.container.css({
      "top": this.top,
      "left": this.left
    });
    this.container.find(".dismiss-notification").click((function(_this) {
      return function() {
        return _this.container.remove();
      };
    })(this));
    $(this.parent).append(this.container);
    this.container.css(css);
  }

  TinyNotification.prototype.dismiss = function() {
    return this.container.find(".dismiss-notification").click();
  };

  return TinyNotification;

})();

//# sourceMappingURL=tiny-notification.js.map
