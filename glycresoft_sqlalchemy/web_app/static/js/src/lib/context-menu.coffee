
contextMenu = (target, options, callback=null) ->
    for item, action of options
        console.log(item, action)
    $(target).off "contextmenu", false
    $(document).off "mousedown", false
    $(target).on "contextmenu", (event) ->
        event.preventDefault()
        handle = $(".context-menu")
        handle.empty()
        if callback?
            callback(handle)
        for item, action of options
            handle.append($("<li></li>").text(item).attr("data-action", item))

        $(".context-menu li").click (e) ->
            handle = $(this)
            console.log(handle)
            action = options[handle.attr("data-action")]
            action(handle)

        $(".context-menu").finish().toggle(100).css(
            {top: event.pageY + 'px', left: event.pageX + 'px'})

$(document).on "mousedown", (e) ->
    if !$(e.target).parents(".context-menu").length > 0
        $(".context-menu").hide(100)


#$ ->
#    $(document).bind("contextmenu", (event) ->
        #
#        # Avoid the real one
#        event.preventDefault();
        #
#        # Show contextmenu
#        $(".context-menu").finish().toggle(100).
        #
#        # In the right position (the mouse)
#        css({
#            top: event.pageY + "px",
#            left: event.pageX + "px"
#        });
#    )
#
#    # If the document is clicked somewhere
#    $(document).bind("mousedown", (e)->
        #
#        # If the clicked element is not the menu
#        if (!$(e.target).parents(".context-menu").length > 0)
            #
#            # Hide it
#            $(".context-menu").hide(100);
#
#    );
#
#
#    # If the menu element is clicked
#    $(".context-menu li").click( (e)->
        #
#        console.log $(this)
#        # This is the triggered action name
#        switch($(this).attr("data-action"))
            #
#            # A case for each action. Your actions here
#            when "first" then alert("first")
#            when "second" then alert("second")
#            when "third" then alert("third")
      #
#        # Hide it AFTER the action was triggered
#        $(".context-menu").hide(100);
#    );
