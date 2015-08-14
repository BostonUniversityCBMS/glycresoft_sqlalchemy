

viewPeakGroupingDatabaseSearchResults = ->
    peptideDetailsModal = undefined

    setup = ->
        $('.protein-match-table tbody tr').click updateProteinChoice
        updateProteinChoice.apply $('.protein-match-table tbody tr')

    updateProteinChoice = ->
        handle = $(this)
        id = handle.attr('data-target')
        $("#chosen-protein-container").html("""<div class="progress"><div class="indeterminate"></div></div>""").fadeIn()
        $.post('/view_database_search_results/protein_composition_view/' + id, GlycReSoft.context).success((doc) ->
            $('#chosen-protein-container').hide()
            $('#chosen-protein-container').html(doc).fadeIn()
            tabs = $('ul.tabs')
            tabs.tabs()
            if GlycReSoft.context['protein-view-active-tab'] != undefined
                console.log GlycReSoft.context['protein-view-active-tab']
                $('ul.tabs').tabs 'select_tab', GlycReSoft.context['protein-view-active-tab']
            else
                $('ul.tabs').tabs 'select_tab', 'protein-overview'
            $('.indicator').addClass 'indigo'
            $('ul.tabs .tab a').click ->
                GlycReSoft.context['protein-view-active-tab'] = $(this).attr('href').slice(1)
            peptideDetailsModal = $('#peptide-detail-modal')
            
            $('.glycopeptide-match-row').click(showGlycopeptideCompositionDetailsModal)

        ).error (error) ->
            console.log arguments

    showGlycopeptideCompositionDetailsModal = ->
        handle = $(this)
        id = handle.attr('data-target')
        console.log(id)
        PartialSource.glycopeptideCompositionDetailsModal {"id": id}, (doc) ->
            peptideDetailsModal.find('.modal-content').html doc
            # Remove any straggler overlays from rapid re-opening of modal
            $(".lean-overlay").remove()
            peptideDetailsModal.openModal()


    unload = ->
        GlycReSoft.removeCurrentLayer()

    setup()
