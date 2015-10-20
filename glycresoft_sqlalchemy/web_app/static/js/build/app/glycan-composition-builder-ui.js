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
