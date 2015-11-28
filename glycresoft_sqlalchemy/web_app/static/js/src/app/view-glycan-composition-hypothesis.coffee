viewGlycanCompositionHypothesis = (hypothesisId) ->
    detailModal = undefined
    displayTable = undefined
    currentPage = 1

    setup = ->
        # $("#save-csv-file").click downloadCSV
        displayTable = $("#composition-table-container") 
        updateCompositionTablePage 1

    setupGlycanCompositionTablePageHandler = (page=1) ->        
        $('.display-table tbody tr').click(->)
        $(':not(.disabled) .next-page').click(-> updateCompositionTablePage(page + 1))
        $(':not(.disabled) .previous-page').click(-> updateCompositionTablePage(page - 1))
        $('.pagination li :not(.active)').click ->
            nextPage = $(@).attr("data-index")
            if nextPage?
                nextPage = parseInt nextPage
                updateCompositionTablePage nextPage

    updateCompositionTablePage = (page=1) ->
        url = "/view_glycan_composition_hypothesis/#{hypothesisId}/#{page}"
        console.log(url)
        GlycReSoft.ajaxWithContext(url).success (doc) ->
            currentPage = page
            displayTable.html doc
            setupGlycanCompositionTablePageHandler page

    setup()
