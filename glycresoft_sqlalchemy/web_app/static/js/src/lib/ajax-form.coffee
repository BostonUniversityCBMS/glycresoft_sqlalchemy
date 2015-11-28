ajaxForm = (formHandle, success, error, transform) ->
    console.log "Ajaxifying ", formHandle
    $(formHandle).on 'submit', (event) ->
        console.log formHandle, "submitting..."
        event.preventDefault()
        handle = $(this)
        if !transform?
            transform = (form) -> new FormData(form)
        url = handle.attr('action')
        method = handle.attr('method')
        data = transform(this)
        encoding = handle.attr('enctype') or 'application/x-www-form-urlencoded; charset=UTF-8'
        ajaxParams = 
            'url': url
            'method': method
            'data': data
            'processData': false
            'contentType': false
            'success': success
            'error': error
        $.ajax ajaxParams


setupAjaxForm = (sourceUrl, container) ->
    container = $(container)
    isModal = container.hasClass('modal')
    $.get(sourceUrl).success (doc) ->
        if isModal
            container.find('.modal-content').html doc
            container.openModal()
            container.find('form').submit (event) ->
                container.closeModal()
    
        else
            container.html doc
    container.find('script').each (i, tag) ->
        tag = $(tag)
        srcURL = tag.attr('src')
        if srcURL != undefined
            $.getScript srcURL
        else
            eval tag.text()
