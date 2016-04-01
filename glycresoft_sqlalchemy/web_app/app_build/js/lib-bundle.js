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

  ActionLayerManager.prototype.setActiveLayerController = function(controller) {
    var layer;
    layer = this.getShowingLayer();
    return layer.setController(controller);
  };

  ActionLayerManager.prototype.getActiveLayerController = function() {
    var layer;
    layer = this.getShowingLayer();
    return layer.controller;
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
    this.controller = null;
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

  ActionLayer.prototype.setController = function(controller) {
    return this.controller = controller;
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

var ColorManager;

ColorManager = (function() {
  function ColorManager(obj) {
    this.store = obj;
  }

  ColorManager.prototype.get = function(name) {
    var b, g, r, rgb, string;
    try {
      rgb = this.store[name];
      r = rgb[0], g = rgb[1], b = rgb[2];
      string = "rgba(" + (r * 255) + "," + (g * 255) + "," + (b * 255) + ",0.5)";
      return string;
    } catch (_error) {
      return this.update;
    }
  };

  ColorManager.prototype.update = function(callback) {
    if (callback == null) {
      callback = function() {
        return {};
      };
    }
    return $.getJSON("/api/colors").success((function(_this) {
      return function(data) {
        _this.store = data;
        return callback(_this);
      };
    })(this));
  };

  return ColorManager;

})();

//# sourceMappingURL=color.js.map

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
      console.log(this, target);
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
  yOffset = 20;
  xOffset = -180;
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

var GlycanComposition;

GlycanComposition = (function() {
  function GlycanComposition(string) {
    this.__order = [];
    this.map = {};
    this.parse(string);
  }

  GlycanComposition.prototype.parse = function(string) {
    var i, len, name, number, part, parts, ref, results;
    parts = string.slice(1, -1).split("; ");
    results = [];
    for (i = 0, len = parts.length; i < len; i++) {
      part = parts[i];
      ref = part.split(":"), name = ref[0], number = ref[1];
      this.__order.push(name);
      results.push(this.map[name] = parseInt(number));
    }
    return results;
  };

  GlycanComposition.prototype.format = function(colorSource) {
    var color, name, number, parts, ref, template;
    parts = [];
    ref = this.map;
    for (name in ref) {
      number = ref[name];
      console.log(name, number);
      if (name === '__order') {
        continue;
      }
      color = colorSource.get(name);
      template = "<span class='monosaccharide-name' style='background-color:" + color + "; padding: 2px;border-radius:2px;'>" + name + " " + number + "</span>";
      parts.push(template);
    }
    return parts.join(' ');
  };

  return GlycanComposition;

})();

//# sourceMappingURL=glycan-composition-parser.js.map

var GlycanComposition;

GlycanComposition = (function() {
  function GlycanComposition(string) {
    this.__order = [];
    this.map = {};
    this.parse(string);
  }

  GlycanComposition.prototype.parse = function(string) {
    var i, len, name, number, part, parts, ref, results;
    parts = string.slice(1, -1).split("; ");
    results = [];
    for (i = 0, len = parts.length; i < len; i++) {
      part = parts[i];
      ref = part.split(":"), name = ref[0], number = ref[1];
      this.__order.push(name);
      results.push(this.map[name] = parseInt(number));
    }
    return results;
  };

  GlycanComposition.prototype.format = function(colorSource) {
    var color, name, number, parts, ref, template;
    parts = [];
    ref = this.map;
    for (name in ref) {
      number = ref[name];
      if (name === '__order') {
        continue;
      }
      color = colorSource.get(name);
      template = "<span class='monosaccharide-name' style='background-color:" + color + "; padding: 2px;border-radius:2px;'>" + name + " " + number + "</span>";
      parts.push(template);
    }
    return parts.join(' ');
  };

  return GlycanComposition;

})();

//# sourceMappingURL=glycan-composition.js.map

var MzIdentMLProteinSelector, getProteinName, getProteinNamesFromMzIdentML, identifyProteomicsFormat;

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
  var id;
  id = /id="([^"]+)"/.exec(line);
  id = id[1];
  return id.split("_").slice(1).join("_");
};

getProteinNamesFromMzIdentML = function(file, callback, nameCallback) {
  var chunksize, fr, offset, proteins, seek;
  fr = new FileReader();
  if (nameCallback == null) {
    nameCallback = function(name) {
      return console.log(name);
    };
  }
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
        if (!proteins[name]) {
          proteins[name] = true;
          nameCallback(name);
        }
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

MzIdentMLProteinSelector = (function() {
  function MzIdentMLProteinSelector(file, listContainer) {
    this.fileObject = file;
    this.container = $(listContainer);
    this.initializeContainer();
  }

  MzIdentMLProteinSelector.prototype.initializeContainer = function() {
    var self, template;
    template = "<div class='display-control'>\n    <a class='toggle-visible-btn right' data-open=\"open\" style='cursor:hand;'>\n        <i class=\"material-icons\">keyboard_arrow_up</i>\n    </a>\n</div>\n<div class='hideable-container'>\n    <div class='row'>\n        <div class='col s4'>\n            <div class='input-field'>\n                <input value='' name=\"protein-regex\" type=\"text\" class=\"validate protein-regex\">\n                <label class=\"active\" for=\"protein-regex\">Protein Pattern</label>\n            </div>\n        </div>\n        <div class='col s2'>\n            <input type='checkbox' id='select-all-proteins-checkbox' name='select-all-proteins-checkbox'/>\n            <label for='select-all-proteins-checkbox'>Select All</label>\n        </div>\n    </div>\n    <div class='row'>\n        <div class='col s8 protein-name-list'>\n\n        </div>\n    </div>\n</div>";
    this.container.html(template);
    this.hideableContainer = this.container.find(".hideable-container");
    this.regex = this.container.find(".protein-regex");
    this.list = this.container.find(".protein-name-list");
    this.toggleVisible = this.container.find(".toggle-visible-btn");
    this.selectAllChecker = this.container.find('#select-all-proteins-checkbox');
    self = this;
    this.toggleVisible.click(function() {
      var handle;
      handle = $(this);
      if (handle.attr("data-open") === "open") {
        self.hideableContainer.hide();
        handle.attr("data-open", "closed");
        return handle.html('<i class="material-icons">keyboard_arrow_down</i>');
      } else if (handle.attr("data-open") === "closed") {
        self.hideableContainer.show();
        handle.attr("data-open", "open");
        return handle.html('<i class="material-icons">keyboard_arrow_up</i>');
      }
    });
    this.regex.change(function(e) {
      var pattern;
      e.preventDefault();
      pattern = $(this).val();
      return self.updateVisibleProteins(pattern);
    });
    this.regex.keydown((function(_this) {
      return function(e) {
        if (e.keyCode === 13) {
          e.preventDefault();
          _this.regex.change();
          return false;
        }
      };
    })(this));
    this.selectAllChecker.click((function(_this) {
      return function(e) {
        if (_this.selectAllChecker.prop("checked")) {
          return _this.container.find("input[type='checkbox']:visible").prop("checked", true);
        } else {
          return _this.container.find("input[type='checkbox']:visible").prop("checked", false);
        }
      };
    })(this));
    return this.load();
  };

  MzIdentMLProteinSelector.prototype.createAddProteinNameToListCallback = function() {
    var callback;
    callback = (function(_this) {
      return function(name) {
        var checker, entryContainer;
        entryContainer = $("<p></p>").css({
          "padding-left": 20,
          "display": 'inline-block',
          "width": 240
        }).addClass('input-field protein-name');
        checker = $("<input />").attr("type", "checkbox").attr("name", name).addClass("protein-name-check");
        entryContainer.append(checker);
        entryContainer.append($("<label></label>").html(name).attr("for", name).click((function() {
          return checker.click();
        })));
        return _this.list.append(entryContainer);
      };
    })(this);
    return callback;
  };

  MzIdentMLProteinSelector.prototype.updateVisibleProteins = function(pattern) {
    var regex;
    regex = new RegExp(pattern, 'i');
    return $('.protein-name').each(function() {
      var handle, name;
      handle = $(this);
      name = handle.find("input").attr("name");
      if (regex.test(name)) {
        return handle.show();
      } else {
        return handle.hide();
      }
    });
  };

  MzIdentMLProteinSelector.prototype.load = function() {
    var callback;
    callback = this.createAddProteinNameToListCallback();
    return getProteinNamesFromMzIdentML(this.fileObject, (function() {}), callback);
  };

  MzIdentMLProteinSelector.prototype.getChosenProteins = function() {
    var a;
    return (function() {
      var i, len, ref, results;
      ref = this.container.find("input.protein-name-check:checked");
      results = [];
      for (i = 0, len = ref.length; i < len; i++) {
        a = ref[i];
        results.push($(a).attr("name"));
      }
      return results;
    }).call(this);
  };

  MzIdentMLProteinSelector.prototype.getAllProteins = function() {
    var a;
    return (function() {
      var i, len, ref, results;
      ref = this.container.find("input.protein-name-check");
      results = [];
      for (i = 0, len = ref.length; i < len; i++) {
        a = ref[i];
        results.push($(a).attr("name"));
      }
      return results;
    }).call(this);
  };

  return MzIdentMLProteinSelector;

})();

//# sourceMappingURL=infer-protein-data-format.js.map


/*
The file name mirroring done by the .file-field group in Materialize is set up on page load.
When these elements are added dynamically, they must be configured manually.

This code is taken from https://github.com/Dogfalo/materialize/blob/master/js/forms.js#L156
 */
var materialCheckbox, materialFileInput, materialRefresh;

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

materialCheckbox = function(selector) {
  return $(selector).click(function(e) {
    var handle, target;
    handle = $(this);
    target = handle.attr("for");
    return $("input[name='" + target + "']").click();
  });
};

//# sourceMappingURL=material-shim.js.map

var PeptideSequence, PeptideSequencePosition, _formatModification;

_formatModification = function(modification, colorSource, long) {
  var color, content;
  if (long == null) {
    long = false;
  }
  color = colorSource.get(modification);
  content = escape(long ? modification : modification[0]);
  return "(<span class='modification-chip' style='background-color:" + color + ";padding-left:1px;padding-right:2px;border-radius:2px;' title='" + (escape(modification)) + "' data-modification=" + (escape(modification)) + ">" + content + "</span>)";
};

PeptideSequencePosition = (function() {
  var formatModification;

  formatModification = function(modification, colorSource, long) {
    var color, content;
    if (long == null) {
      long = false;
    }
    color = colorSource.get(modification);
    content = long ? modification : modification[0];
    return "(<span class='modification-chip' style='background-color:" + color + ";padding-left:1px;padding-right:2px;border-radius:2px;' title='" + modification + "' data-modification=" + modification + ">" + content + "</span>)";
  };

  function PeptideSequencePosition(residue, modifications) {
    if (modifications == null) {
      modifications = [];
    }
    this.residue = residue;
    this.modifications = modifications;
  }

  PeptideSequencePosition.prototype.format = function(colorSource) {
    var modifications;
    modifications = this.modifications.map(function(modification) {
      return formatModification(modification, colorSource);
    }).join('');
    return this.residue + modifications;
  };

  return PeptideSequencePosition;

})();

PeptideSequence = (function() {
  var parser;

  parser = function(sequence) {
    var cTerm, chunks, currentAA, currentMod, currentMods, glycan, i, m, mods, nTerm, nextChr, p, parenLevel, state;
    state = "start";
    nTerm = nTerm || "H";
    cTerm = cTerm || "OH";
    mods = [];
    chunks = [];
    glycan = "";
    currentAA = "";
    currentMod = "";
    currentMods = [];
    parenLevel = 0;
    i = 0;
    while (i < sequence.length) {
      nextChr = sequence[i];
      if (nextChr === "(") {
        if (state === "aa") {
          state = "mod";
          parenLevel += 1;
        } else if (state === "start") {
          state = "nTerm";
          parenLevel += 1;
        } else {
          parenLevel += 1;
          if (!((state === "nTerm" || state === "cTerm") && parenLevel === 1)) {
            currentMod += nextChr;
          }
        }
      } else if (nextChr === ")") {
        if (state === "aa") {
          throw new Exception("Invalid Sequence. ) found outside of modification, Position {0}. {1}".format(i, sequence));
        } else {
          parenLevel -= 1;
          if (parenLevel === 0) {
            mods.push(currentMod);
            currentMods.push(currentMod);
            if (state === "mod") {
              state = 'aa';
              if (currentAA === "") {
                chunks.slice(-1)[1] = chunks.slice(-1)[1].concat(currentMods);
              } else {
                chunks.push([currentAA, currentMods]);
              }
            } else if (state === "nTerm") {
              if (sequence[i + 1] !== "-") {
                throw new Exception("Malformed N-terminus for " + sequence);
              }
              nTerm = currentMod;
              state = "aa";
              i += 1;
            } else if (state === "cTerm") {
              cTerm = currentMod;
            }
            currentMods = [];
            currentMod = "";
            currentAA = "";
          } else {
            currentMod += nextChr;
          }
        }
      } else if (nextChr === "|") {
        if (state === "aa") {
          throw new Exception("Invalid Sequence. | found outside of modification");
        } else {
          currentMods.push(currentMod);
          mods.push(currentMod);
          currentMod = "";
        }
      } else if (nextChr === "{") {
        if (state === 'aa' || (state === "cTerm" && parenLevel === 0)) {
          glycan = sequence.slice(i);
          break;
        } else {
          currentMod += nextChr;
        }
      } else if (nextChr === "-") {
        if (state === "aa") {
          state = "cTerm";
          if (currentAA !== "") {
            currentMods.push(currentMod);
            chunks.push([currentAA, currentMods]);
            currentMod = "";
            currentMods = [];
            currentAA = "";
          }
        } else {
          currentMod += nextChr;
        }
      } else if (state === "start") {
        state = "aa";
        currentAA = nextChr;
      } else if (state === "aa") {
        if (currentAA !== "") {
          currentMods.push(currentMod);
          chunks.push([currentAA, currentMods]);
          currentMod = "";
          currentMods = [];
          currentAA = "";
        }
        currentAA = nextChr;
      } else if (state === "nTerm" || state === "mod" || state === "cTerm") {
        currentMod += nextChr;
      } else {
        throw new Exception("Unknown Tokenizer State", currentAA, currentMod, i, nextChr);
      }
      i += 1;
    }
    if (currentAA !== "") {
      chunks.push([currentAA, currentMod]);
    }
    if (currentMod !== "") {
      mods.push(currentMod);
    }
    if (glycan !== "") {
      glycan = new GlycanComposition(glycan);
    } else {
      glycan = null;
    }
    chunks = (function() {
      var j, len, results;
      results = [];
      for (j = 0, len = chunks.length; j < len; j++) {
        p = chunks[j];
        results.push(new PeptideSequencePosition(p[0], (function() {
          var k, len1, ref, results1;
          ref = p[1];
          results1 = [];
          for (k = 0, len1 = ref.length; k < len1; k++) {
            m = ref[k];
            if (m !== "") {
              results1.push(m);
            }
          }
          return results1;
        })()));
      }
      return results;
    })();
    return [chunks, mods, glycan, nTerm, cTerm];
  };

  function PeptideSequence(string) {
    var mods, ref;
    ref = parser(string), this.sequence = ref[0], mods = ref[1], this.glycan = ref[2], this.nTerm = ref[3], this.cTerm = ref[4];
  }

  PeptideSequence.prototype.format = function(colorSource, includeGlycan) {
    var cTerm, glycan, nTerm, position, positions, sequence;
    if (includeGlycan == null) {
      includeGlycan = true;
    }
    positions = [];
    if (this.nTerm === "H") {
      nTerm = "";
    } else {
      nTerm = _formatModification(this.nTerm, colorSource).slice(1, -1) + '-';
      positions.push(nTerm);
    }
    positions = positions.concat((function() {
      var j, len, ref, results;
      ref = this.sequence;
      results = [];
      for (j = 0, len = ref.length; j < len; j++) {
        position = ref[j];
        results.push(position.format(colorSource));
      }
      return results;
    }).call(this));
    if (this.cTerm === "OH") {
      cTerm = "";
    } else {
      cTerm = '-' + (_formatModification(this.cTerm, colorSource).slice(1, -1));
      positions.push(cTerm);
    }
    sequence = positions.join("");
    if (includeGlycan) {
      glycan = this.glycan.format(colorSource);
      sequence += ' ' + glycan;
    }
    return sequence;
  };

  return PeptideSequence;

})();

//# sourceMappingURL=peptide-sequence.js.map

var TinyNotification, tinyNotify;

TinyNotification = (function() {
  TinyNotification.prototype.template = "<div class='notification-container'>\n    <div class='clearfix dismiss-container'>\n        <a class='dismiss-notification mdi-content-clear'></a>\n    </div>\n    <div class='notification-content'>\n    </div>\n</div>";

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

tinyNotify = function(top, left, message, parent, css) {
  if (parent == null) {
    parent = 'body';
  }
  if (css == null) {
    css = {};
  }
  return new TinyNotification(top, left, message, parent, css);
};

//# sourceMappingURL=tiny-notification.js.map
