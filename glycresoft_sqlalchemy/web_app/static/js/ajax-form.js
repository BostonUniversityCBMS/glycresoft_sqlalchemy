function ajaxForm(formHandle, success, error){
    $(formHandle).on("submit", function(event){
        event.preventDefault();
        var handle = $(this)
        var url = handle.attr("action")
        var method = handle.attr("method")
        var data = new FormData(this)
        var encoding = handle.attr('enctype') || "application/x-www-form-urlencoded; charset=UTF-8"
        var ajaxParams = {
            "url": url,
            "method": method,
            "data": data,
            "processData": false,
            "contentType": false,
            "success": success,
            "error": error
        }
        $.ajax(ajaxParams)
    })
}


function setupAjaxForm(sourceUrl, container){
    container = $(container)
    var isModal = container.hasClass('modal');
    $.get(sourceUrl).success(function(doc){
        if(isModal){
            container.find(".modal-content").html(doc)
            container.openModal();
            container.find("form").submit(function(event){
                container.closeModal();
            })
        } else {
            container.html(doc)
        }
    })
    container.find("script").each(function(i, tag){
        tag = $(tag)
        var srcURL = tag.attr("src")
        if(srcURL !== undefined){
            $.getScript(srcURL)
        } else {
            eval(tag.text())
        }
    })
}


