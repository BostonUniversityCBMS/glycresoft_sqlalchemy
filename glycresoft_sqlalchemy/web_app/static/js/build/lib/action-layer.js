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

  ActionLayerManager.prototype.addLayer = function(options) {
    return new ActionLayer(this, options);
  };

  ActionLayerManager.prototype.removeLayer = function(id) {
    this.layers[id].container.remove();
    return delete this.layers[id];
  };

  return ActionLayerManager;

})(EventEmitter);

ActionLayer = (function() {
  function ActionLayer(manager, options) {
    this.manager = manager;
    this.options = options;
    this.contentURL = options.contentURL;
    if (!options.container) {
      this.id = options.name || 'action-layer-' + manager.incLayerCounter();
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
    self = this;
    return $.get(this.contentURL).success(function(doc) {
      self.container.hide();
      self.container.html(doc);
      return self.container.find('script').each(function(i, tag) {
        var srcURL;
        tag = $(tag);
        srcURL = tag.attr('src');
        if (srcURL !== void 0) {
          return $.getScript(srcURL);
        } else {
          return eval(tag.text());
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

  return ActionLayer;

})();

//# sourceMappingURL=action-layer.js.map
