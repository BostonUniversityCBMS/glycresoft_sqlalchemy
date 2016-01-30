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
