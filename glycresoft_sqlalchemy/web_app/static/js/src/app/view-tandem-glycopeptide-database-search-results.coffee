doZoom = ->
    svg = d3.select("svg g")
    zoom = ()->
        svg.attr("transform", "translate(#{d3.event.translate})scale(#{d3.event.scale})")
    d3.select("svg g").call(d3.behavior.zoom().scaleExtent([1, 8]).on("zoom", zoom))


viewTandemGlycopeptideDatabaseSearchResults = ->
    peptideDetailsModal = undefined

    setup = ->
        $('.protein-match-table tbody tr').click updateProteinChoice

        last_id = GlycReSoft.context['protein_id']
        last_selector = '''.protein-match-table tbody tr[data-target="''' + last_id + '''"]'''
        console.log(last_selector)
        handle = $(last_selector)
        console.log(handle)
        if handle.length != 0
            updateProteinChoice.apply handle
        else    
            updateProteinChoice.apply $('.protein-match-table tbody tr')
        $(".tooltipped").tooltip()
        $("#save-csv-file").click downloadCSV

    initGlycopeptideOverviewPlot = ->
        $('svg .glycopeptide').customTooltip glycopeptideTooltipCallback, 'protein-view-tooltip'
        $('svg .modification path').customTooltip modificationTooltipCallback, 'protein-view-tooltip'


    glycopeptideTooltipCallback = (handle) ->
        template = '<div>
        <span><b>MS2 Score:</b> {ms2-score}</span><br>
        <span><b>q-value:</b> {q-value}</span><br>
        <b>{sequence}</b>
        </div>'
        template.format
            'sequence': handle.attr('data-sequence')
            'ms2-score': handle.attr('data-ms2-score')
            'q-value': handle.attr('data-q-value')

    modificationTooltipCallback = (handle) ->
        template = '
        <div>
        <span>{value}</span>
        </div>'
        value = handle.parent().attr('data-modification-type')
        if value == 'HexNAc'
            sequence = $('#' + handle.parent().attr('data-parent')).attr('data-sequence')
            value = 'HexNAc - Glycosylation: ' + sequence.split(/(\[|\{)/).slice(1).join('')
        template.format 'value': value

    updateProteinChoice = ->
        handle = $(this)
        id = handle.attr('data-target')
        $("#chosen-protein-container").html("""<div class="progress"><div class="indeterminate"></div></div>""").fadeIn()
        $.ajax '/view_database_search_results/protein_view/' + id,
                data: JSON.stringify({"context": GlycReSoft.context, "settings": GlycReSoft.settings}),
                contentType: "application/json"
                type: 'POST'
                success: (doc) ->
                    $('#chosen-protein-container').hide()
                    $('#chosen-protein-container').html(doc).fadeIn()
                    initGlycopeptideOverviewPlot()
                    tabs = $('ul.tabs')
                    tabs.tabs()
                    if GlycReSoft.context['protein-view-active-tab'] != undefined
                        console.log GlycReSoft.context['protein-view-active-tab']
                        $('ul.tabs').tabs 'select_tab', GlycReSoft.context['protein-view-active-tab']
                    else
                        $('ul.tabs').tabs 'select_tab', 'protein-overview'
                    $('ul.tabs .tab a').click ->
                        GlycReSoft.context['protein-view-active-tab'] = $(this).attr('href').slice(1)
                    $('.indicator').addClass 'indigo'
                    $('.glycopeptide-match-row').click showGlycopeptideDetailsModal
                    peptideDetailsModal = $('#peptide-detail-modal')
                    GlycReSoft.context['protein_id'] = id
                error: (error) ->
                    console.log arguments

    getGlycopeptideMatchDetails = (id, callback) ->
        $.get '/api/glycopeptide_match/' + id, callback

    showGlycopeptideDetailsModal = ->
        handle = $(this)
        id = handle.attr('data-target')
        $.get('/view_database_search_results/view_glycopeptide_details/' + id).success (doc) ->
            peptideDetailsModal.find('.modal-content').html doc
            # Remove any straggler overlays from rapid re-opening of modal
            $(".lean-overlay").remove()
            peptideDetailsModal.openModal()

    downloadCSV = ->
        handle = $(this)
        id = handle.attr('data-target')
        $.ajax "/view_database_search_results/export_csv/" + id,
                data: JSON.stringify({"context": GlycReSoft.context, "settings": GlycReSoft.settings}),
                contentType: "application/json"
                type: 'POST'

    setup()
