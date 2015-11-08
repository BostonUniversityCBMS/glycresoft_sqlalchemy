

viewPeakGroupingDatabaseSearchResults = ->
    glycopeptideDetailsModal = undefined
    glycopeptideTable = undefined
    currentPage = 1
    currentProtein = undefined

    setup = ->
        $('.protein-match-table tbody tr').click updateProteinChoice
        updateProteinChoice.apply $('.protein-match-table tbody tr')
        console.log "glycopeptideTable", glycopeptideTable

    setupGlycopeptideCompositionTablePageHandlers = (page=1) ->
        $('.glycopeptide-match-row').click(showGlycopeptideCompositionDetailsModal)
        $(':not(.disabled) .next-page').click(-> updateGlycopeptideCompositionTablePage(page + 1))
        $(':not(.disabled) .previous-page').click(-> updateGlycopeptideCompositionTablePage(page - 1))
        $('.pagination li :not(.active)').click ->
            nextPage = $(@).attr("data-index")
            if nextPage?
                nextPage = parseInt nextPage
                updateGlycopeptideCompositionTablePage nextPage

    updateGlycopeptideCompositionTablePage = (page=1) ->
        url = "/view_database_search_results/glycopeptide_matches_composition_table/#{currentProtein}/#{page}"
        console.log(url)
        GlycReSoft.ajaxWithContext(url).success (doc) ->
            currentPage = page
            glycopeptideTable.html doc
            setupGlycopeptideCompositionTablePageHandlers page

    updateProteinChoice = ->
        handle = $(this)
        currentProtein = id = handle.attr('data-target')
        console.log glycopeptideDetailsModal
        $("#chosen-protein-container").html("""<div class="progress"><div class="indeterminate"></div></div>""").fadeIn()
        $.post('/view_database_search_results/protein_composition_view/' + id, GlycReSoft.context).success((doc) ->
            $('#chosen-protein-container').hide()
            $('#chosen-protein-container').html(doc).fadeIn()
            tabs = $('ul.tabs')
            tabs.tabs()
            GlycReSoft.context["current_protein"] = id
            if GlycReSoft.context['protein-view-active-tab'] != undefined
                console.log GlycReSoft.context['protein-view-active-tab']
                $('ul.tabs').tabs 'select_tab', GlycReSoft.context['protein-view-active-tab']
            else
                $('ul.tabs').tabs 'select_tab', 'protein-overview'
            $('.indicator').addClass 'indigo'
            $('ul.tabs .tab a').click ->
                GlycReSoft.context['protein-view-active-tab'] = $(this).attr('href').slice(1)
            glycopeptideDetailsModal = $('#peptide-detail-modal')
            glycopeptideTable = $("#glycopeptide-table")
            setupGlycopeptideCompositionTablePageHandlers(1)
        ).error (error) ->
            console.log arguments

    showGlycopeptideCompositionDetailsModal = ->
        handle = $(this)
        id = handle.attr('data-target')
        console.log(glycopeptideDetailsModal)
        console.log(id)
        PartialSource.glycopeptideCompositionDetailsModal {"id": id}, (doc) ->
            glycopeptideDetailsModal.find('.modal-content').html doc
            $(".lean-overlay").remove()
            glycopeptideDetailsModal.openModal()

    unload = ->
        GlycReSoft.removeCurrentLayer()

    setup()
