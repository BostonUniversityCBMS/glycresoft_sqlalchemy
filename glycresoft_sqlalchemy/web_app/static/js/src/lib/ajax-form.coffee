ajaxForm = (formHandle, success, error) ->
    $(formHandle).on 'submit', (event) ->
        event.preventDefault()
        handle = $(this)
        url = handle.attr('action')
        method = handle.attr('method')
        data = new FormData(this)
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
