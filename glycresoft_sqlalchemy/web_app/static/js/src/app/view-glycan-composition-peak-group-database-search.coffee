doZoom = (selector) ->
        svg = d3.select(selector)
        zoom = ()->
            svg.attr("transform", "translate(#{d3.event.translate})scale(#{d3.event.scale})")
        d3.select(selector).call(d3.behavior.zoom().scaleExtent([1, 8]).on("zoom", zoom))

viewGlycanCompositionPeakGroupingDatabaseSearchResults = ->
    glycanDetailsModal = undefined
    glycanTable = undefined
    currentPage = 1

    setup = ->
        updateView()
        $("#save-csv-file").click downloadCSV

    setupGlycanCompositionTablePageHandlers = (page=1) ->
        $('.glycan-match-row').click(showGlycanCompositionDetailsModal)
        $(':not(.disabled) .next-page').click(-> updateGlycanCompositionTablePage(page + 1))
        $(':not(.disabled) .previous-page').click(-> updateGlycanCompositionTablePage(page - 1))
        $('.pagination li :not(.active)').click ->
            nextPage = $(@).attr("data-index")
            if nextPage?
                nextPage = parseInt nextPage
                updateGlycanCompositionTablePage nextPage

    updateGlycanCompositionTablePage = (page=1) ->
        url = "/view_database_search_results/glycan_composition_match_table/#{page}"
        console.log(url)
        GlycReSoft.ajaxWithContext(url).success (doc) ->
            currentPage = page
            glycanTable.html doc
            setupGlycanCompositionTablePageHandlers page

    updateView = ->
        handle = $(this)

        $("#content-container").html("""<div class="progress"><div class="indeterminate"></div></div>""").fadeIn()
        GlycReSoft.ajaxWithContext('/view_database_search_results/results_view/').success((doc) ->
            $('#content-container').hide()
            $('#content-container').html(doc).fadeIn()
            tabs = $('ul.tabs')
            tabs.tabs()
            if GlycReSoft.context['view-active-tab'] != undefined
                console.log GlycReSoft.context['view-active-tab']
                $('ul.tabs').tabs 'select_tab', GlycReSoft.context['view-active-tab']
            else
                $('ul.tabs').tabs 'select_tab', 'glycome-overview'
            $('.indicator').addClass 'indigo'
            $('ul.tabs .tab a').click ->
                GlycReSoft.context['view-active-tab'] = $(this).attr('href').slice(1)
            glycanDetailsModal = $ '#glycan-detail-modal'
            glycanTable = $ "#glycan-table"
            # doZoom("#glycome-overview svg g")
            updateGlycanCompositionTablePage 1
        ).error (error) ->
            console.log arguments

    showGlycanCompositionDetailsModal = ->
        handle = $(this)
        id = handle.attr('data-target')
        console.log id
        PartialSource.glycanCompositionDetailsModal {"id": id}, (doc) ->
            glycanDetailsModal.find('.modal-content').html doc
            $(".lean-overlay").remove()
            glycanDetailsModal.openModal()

    unload = ->
        GlycReSoft.removeCurrentLayer()


    downloadCSV = ->
        handle = $(this)
        id = handle.attr('data-target')
        $.ajax "/view_database_search_results/export_csv/" + id,
                data: JSON.stringify({"context": GlycReSoft.context, "settings": GlycReSoft.settings}),
                contentType: "application/json"
                type: 'POST'

    setup()
