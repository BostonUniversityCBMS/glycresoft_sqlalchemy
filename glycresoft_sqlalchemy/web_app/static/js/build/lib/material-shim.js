
/*
The file name mirroring done by the .file-field group in Materialize is set up on page load.
When these elements are added dynamically, they must be configured manually.

This code is taken from https://github.com/Dogfalo/materialize/blob/master/js/forms.js#L156
 */
var materialFileInput, materialRefresh;

materialRefresh = function() {
  materialFileInput();
  Materialize.updateTextFields();
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
