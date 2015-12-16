
contextMenu = (target, options, callback=null) ->
    # Clear existing menu configuration
    $(target).off "contextmenu", false
    $(document).off "mousedown", false

    # Attach new contextmenu instructions
    $(target).on "contextmenu", (event) ->
        event.preventDefault()
        handle = $(".context-menu")
        handle.empty()
        # If a callback is provided to configure the menu,
        # execute it.
        if callback?
            callback(handle)
        # For each item in `options`, create an
        # entry in the menu and bind its paired action
        for item, action of options
            handle.append($("<li></li>").text(item).attr("data-action", item))

        # Set up event handlers for each context menu item
        # to invoke clicked items' actions on the owner of
        # the menu
        $(".context-menu li").click (e) ->
            handle = $(this)
            console.log this, target
            action = options[handle.attr("data-action")]
            action.apply target

        # Finally display the context menu at the site of
        # the event.
        $(".context-menu").finish().toggle(100).css(
            {top: event.pageY + 'px', left: event.pageX + 'px'})

# If there is a mousedown event not inside of a context menu,
# hide any extant menus.
$(document).on "mousedown", (e) ->
    if !$(e.target).parents(".context-menu").length > 0
        $(".context-menu").hide(100)
