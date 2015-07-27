$(function(){
    var $body = $("body")
    var $tooltip = $("<div></div>").hide().css({"position": "absolute", "z-index": "10"})
    $body.append($tooltip)


    function openTooltip(event){
        var handle = $(this)
        var content = handle.data("tooltip-content")
        if (typeof content === "function"){
            content = content(handle)
        }
        content = content === undefined ? "This is a simple tooltip" : content
        $tooltip.html(content)
        .addClass(handle.data("tooltip-css-class")).css(
            "top", (event.pageY-10)+"px").css("left",(event.pageX+10)+"px")
        .show()
    }

    function closeTooltip(event){
        var handle = $(this)
        $tooltip.html("")
        .removeClass(handle.data("tooltip-css-class"))
        .hide()
    }

    jQuery.fn.customTooltip = function(content, cssClass){
        var handle = $(this)
        if (content !== undefined){
            handle.data("tooltip-content", content)
        }
        if (cssClass !== undefined){
            handle.data("tooltip-css-class", cssClass)
        }
        handle.hover(openTooltip, closeTooltip)
    }
})
