$ ->
    yOffset = 20 # -3
    xOffset = -180
    $body = $('body')
    $tooltip = $('<div></div>').hide().css(
        'position': 'absolute'
        'z-index': '10')

    openTooltip = (event) ->
        handle = $(this)
        content = handle.data('tooltip-content')
        if typeof content == 'function'
            content = content(handle)
        content = if content == undefined then 'This is a simple tooltip' else content
        $tooltip.html(content).addClass(
            handle.data('tooltip-css-class')).css(
                'top', event.pageY + yOffset + 'px').css(
                'left', event.pageX + xOffset + 'px').show()
        return

    closeTooltip = (event) ->
        handle = $(this)
        $tooltip.html('').removeClass(handle.data('tooltip-css-class')).hide()
        return

    $body.append $tooltip

    jQuery.fn.customTooltip = (content, cssClass) ->
        handle = $(this)
        if content != undefined
            handle.data 'tooltip-content', content
        if cssClass != undefined
            handle.data 'tooltip-css-class', cssClass
        handle.hover openTooltip, closeTooltip
        return

    return
