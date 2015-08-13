###
The file name mirroring done by the .file-field group in Materialize is set up on page load.
When these elements are added dynamically, they must be configured manually.

This code is taken from https://github.com/Dogfalo/materialize/blob/master/js/forms.js#L156
###

materialRefresh = ->
    try
        $('select').material_select();
    try
        materialFileInput()
    try
        Materialize.updateTextFields()
    return

materialFileInput = ->
    $(document).on 'change', '.file-field input[type="file"]', ->
        file_field = $(this).closest('.file-field')
        path_input = file_field.find('input.file-path')
        files = $(this)[0].files
        file_names = []
        i = 0
        while i < files.length
            file_names.push files[i].name
            i++
        path_input.val file_names.join(', ')
        path_input.trigger 'change'
        return
    return
