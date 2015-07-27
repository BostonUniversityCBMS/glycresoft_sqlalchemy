class ViewDatabaseSearchResults extends EventEmitter
    constructor: ->
        $('.protein-match-table tbody tr').click updateProteinChoice
        updateProteinChoice.apply $('.protein-match-table tbody tr'), null
        return

    initGlycopeptideTooltip: ->
        $('svg .glycopeptide').customTooltip @glycopeptideTooltipCallback, 'protein-view-tooltip'
        $('svg .modification path').customTooltip @modificationTooltipCallback, 'protein-view-tooltip'
        return

    glycopeptideTooltipCallback: (handle) ->
        template = '<div>
        <span><b>MS2 Score:</b> {ms2-score}</span><br>
        <span><b>q-value:</b> {q-value}</span><br>
        <b>{sequence}</b>
        </div>'
        template.format
            'sequence': handle.attr('data-sequence')
            'ms2-score': handle.attr('data-ms2-score')
            'q-value': handle.attr('data-q-value')

    modificationTooltipCallback: (handle) ->
        template = '
        <div>
        <span>{value}</span>
        </div>'
        value = handle.parent().attr('data-modification-type')
        if value == 'HexNAc'
            sequence = $('#' + handle.parent().attr('data-parent')).attr('data-sequence')
            value = 'HexNAc - Glycosylation: ' + sequence.split(/(\[|\{)/).slice(1).join('')
        template.format 'value': value

peptideDetailsModal = undefined

updateProteinChoice = ->
    handle = $(this)
    id = handle.attr('data-target')
    console.log id
    $('#chosen-protein-container').fadeOut()
    $.get('/view_database_search_results/protein_view/' + id).success((doc) ->
        $('#chosen-protein-container').html(doc).fadeIn()
        initGlycopeptideTooltip()
        tabs = $('ul.tabs')
        tabs.tabs()
        if GlycReSoft.context['protein-view-active-tab'] != undefined
            console.log GlycReSoft.context['protein-view-active-tab']
            $('ul.tabs').tabs 'select_tab', GlycReSoft.context['protein-view-active-tab']
        else
            $('ul.tabs').tabs 'select_tab', 'protein-overview'
        $('ul.tabs .tab a').click ->
            GlycReSoft.context['protein-view-active-tab'] = $(this).attr('href').slice(1)
            return
        $('.indicator').addClass 'indigo'
        $('.glycopeptide-match-row').click showGlycopeptideDetailsModal
        peptideDetailsModal = $('#peptide-detail-modal')
    ).error (error) ->
        console.log arguments

getGlycopeptideMatchDetails = (id, callback) ->
    $.get '/api/glycopeptide_match/' + id, callback
    return

showGlycopeptideDetailsModal = ->
    handle = $(this)
    id = handle.attr('data-target')
    $.get('/view_database_search_results/view_glycopeptide_details/' + id).success (doc) ->
        peptideDetailsModal.find('.modal-content').html doc
        peptideDetailsModal.openModal()
