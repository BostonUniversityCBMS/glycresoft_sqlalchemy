var MassShiftInputWidget;

MassShiftInputWidget = (function() {
  var addEmptyRowOnEdit, counter, template;
  template = "<div class='mass-shift-row row'>\n    <div class='input-field col s3'>\n        <label for='mass_shift_name'>Mass Shift Name</label>\n        <input class='mass-shift-name' type='text' name='mass_shift_name' placeholder='Name'>\n    </div>\n    <div class='input-field col s3'>\n        <label for='mass_shift_mass_delta'>Mass &Delta;</label>\n        <input class='mass-delta' type='number' name='mass_shift_mass_delta' step=\"0.0001\" placeholder='Mass Shift'>\n    </div>\n    <div class='input-field col s3'>\n        <label for='mass_shift_max_count'>Maximum Count</label>    \n        <input class='max-count' type='number' min='0' placeholder='Maximum Count' name='mass_shift_max_count'>\n    </div>\n</div>";
  counter = 0;
  addEmptyRowOnEdit = function(container, addHeader) {
    var callback, row;
    if (addHeader == null) {
      addHeader = true;
    }
    container = $(container);
    if (addHeader) {
      row = $(template);
    } else {
      row = $(template);
      row.find("label").remove();
    }
    container.append(row);
    row.data("counter", ++counter);
    callback = function(event) {
      if (row.data("counter") === counter) {
        addEmptyRowOnEdit(container, false);
      }
      return $(this).parent().find("label").removeClass("active");
    };
    return row.find("input").change(callback);
  };
  return addEmptyRowOnEdit;
})();

//# sourceMappingURL=mass-shift-ui.js.map
