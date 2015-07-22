function initGlycopeptideTooltip(){
    $("svg .glycopeptide").map(function(i, elem){
        elem = $(elem)
        elem.attr("data-tooltip", elem.attr("data-sequence"))
    })
    $("svg .glycopeptide").tooltip()

/*    $("svg .modification").map(function(i, elem){
        elem = $(elem)
        elem.find("path").attr("data-tooltip", elem.attr("data-modification-type"))        
    })
    $("svg .modification path").tooltip()
*/}

$(function(){
    $("select.protein-selection").change(function(){
        var handle = $(this)
        console.log(handle.val())
        $.get("/view_database_search_results/protein_view/" + handle.val()).success(
            function(doc){
                $("#chosen-protein-container").html(doc)
                initGlycopeptideTooltip()
            }).error(function(error){
                console.log(arguments)
            })
    })
})

