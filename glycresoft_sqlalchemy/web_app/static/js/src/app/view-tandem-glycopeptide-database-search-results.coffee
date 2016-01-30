TandemGlycopeptideDatabaseSearchResultsController = (->
    glycopeptideTooltipCallback = (handle) ->
        template = '''<div><table>
        <tr><td style='padding:3px;'><b>MS2 Score:</b> {ms2-score}</td><td style='padding:3px;'><b>Mass:</b> {calculated-mass}</td></tr>
        <tr><td style='padding:3px;'><b>q-value:</b> {q-value}</td><td style='padding:3px;'><b>Spectrum Matches:</b> {spectra-count}</td></tr>
        </table>
        <span>{sequence}</span>
        </div>'''

        template.format
            'sequence': new PeptideSequence(handle.attr('data-sequence')).format(GlycReSoft.colors)
            'ms2-score': parseFloat(handle.attr('data-ms2-score')).toFixed(4)
            'q-value': handle.attr('data-q-value')
            "calculated-mass": parseFloat(handle.attr("data-calculated-mass")).toFixed(4)
            "spectra-count": handle.attr("data-spectra-count")

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

    getGlycopeptideMatchDetails = (id, callback) ->
        $.get '/api/glycopeptide_match/' + id, callback

    downloadCSV = ->
        handle = $(this)
        id = handle.attr('data-target')
        $.ajax "/view_database_search_results/export_csv/" + id,
                data: JSON.stringify({"context": GlycReSoft.context, "settings": GlycReSoft.settings}),
                contentType: "application/json"
                type: 'POST'

    class TandemGlycopeptideDatabaseSearchResultsController
        constructor: ->
            @currentPage = 1
            @peptideDetailsModal = undefined
            @glycopeptideTable = undefined
            @currentProtein = undefined

            @setup()


        setup: ->
            updateProteinChoice = @updateProteinChoiceCallback()
            $('.protein-match-table tbody tr').click updateProteinChoice

            last_id = GlycReSoft.context['protein_id']
            last_selector = '''.protein-match-table tbody tr[data-target="''' + last_id + '''"]'''
            handle = $(last_selector)

            if handle.length != 0
                updateProteinChoice.apply handle
            else    
                updateProteinChoice.apply $($('.protein-match-table tbody tr')[0])
            $(".tooltipped").tooltip()
            $("#save-csv-file").click downloadCSV

        setupGlycopeptideTablePageHandlers: (page=1) ->
            self = @
            $('.glycopeptide-match-row').click ->
                textSelection = window.getSelection()
                if not textSelection.toString()
                    self.showGlycopeptideDetailsModalCallback().apply @
            $(':not(.disabled) .next-page').click(-> self.updateGlycopeptideTablePage(page + 1))
            $(':not(.disabled) .previous-page').click(-> self.updateGlycopeptideTablePage(page - 1))
            $('.pagination li :not(.active)').click ->
                nextPage = $(@).attr("data-index")
                if nextPage?
                    nextPage = parseInt nextPage
                    self.updateGlycopeptideTablePage nextPage

        updateGlycopeptideTablePage: (page=1) ->
            url = "/view_database_search_results/glycopeptide_match_table/#{@currentProtein}/#{page}"
            GlycReSoft.ajaxWithContext(url).success (doc) =>
                @currentPage = page
                @glycopeptideTable.html doc
                @setupGlycopeptideTablePageHandlers page

        initGlycopeptideOverviewPlot: ->
            glycopeptide = $('svg .glycopeptide')
            glycopeptide.customTooltip glycopeptideTooltipCallback, 'protein-view-tooltip'
            self = @
            glycopeptide.hover(
                (event) ->
                    handle = $ @
                    baseColor = handle.find("path").css("fill")
                    newColor = '#74DEC5'
                    handle.data("baseColor", baseColor)
                    handle.find("path").css("fill", newColor)
                (event) ->
                    handle = $ @
                    handle.find("path").css("fill", handle.data("baseColor"))
                )
            glycopeptide.click (event) ->
                handle = $ @
                id = handle.data("record-id")
                $.get('/view_database_search_results/view_glycopeptide_details/' + id).success (doc) ->
                    self.peptideDetailsModal.find('.modal-content').html doc
                    # Remove any straggler overlays from rapid re-opening of modal
                    $(".lean-overlay").remove()
                    self.peptideDetailsModal.openModal()
            $('svg .modification path').customTooltip modificationTooltipCallback, 'protein-view-tooltip'

        updateProteinChoiceCallback: ->
            self = this
            callback = ->
                handle = $(this)
                $('.active-row').removeClass("active-row")
                handle.addClass("active-row")
                id = handle.attr('data-target')
                self.currentProtein = id
                $('#chosen-protein-container').fadeOut()
                $("#loading-top-level-chosen-protein-container").fadeIn()
                $.ajax '/view_database_search_results/protein_view/' + id,
                        data: JSON.stringify({"context": GlycReSoft.context, "settings": GlycReSoft.settings}),
                        contentType: "application/json"
                        type: 'POST'
                        success: (doc) ->
                            $('#chosen-protein-container').hide()
                            $("#loading-top-level-chosen-protein-container").fadeOut()
                            $('#chosen-protein-container').html(doc).fadeIn()
                            tabs = $('ul.tabs')
                            GlycReSoft.ajaxWithContext(
                                "/view_database_search_results/protein_view/#{id}/protein_overview_panel").success (svg) ->
                                    $("#protein-overview").html(svg)
                                    self.initGlycopeptideOverviewPlot()
                            GlycReSoft.ajaxWithContext(
                                "/view_database_search_results/protein_view/#{id}/microheterogeneity_plot_panel").success (svgGal) ->
                                    $("#site-distribution").html(svgGal)

                            tabs.tabs()
                            if GlycReSoft.context['protein-view-active-tab'] != undefined
                                $('ul.tabs').tabs 'select_tab', GlycReSoft.context['protein-view-active-tab']
                            else
                                $('ul.tabs').tabs 'select_tab', 'protein-overview'
                            $('ul.tabs .tab a').click ->
                                GlycReSoft.context['protein-view-active-tab'] = $(this).attr('href').slice(1)
                            $('.indicator').addClass 'indigo'
                            self.glycopeptideTable = $ "#glycopeptide-table"
                            self.updateGlycopeptideTablePage 1
                            self.peptideDetailsModal = $('#peptide-detail-modal')
                            GlycReSoft.context['protein_id'] = id
                        error: (error) ->
                            console.log arguments

        showGlycopeptideDetailsModalCallback: ->
            self = @
            callback = ->
                handle = $(this)
                id = handle.attr('data-target')
                $.get('/view_database_search_results/view_glycopeptide_details/' + id).success (doc) ->
                    self.peptideDetailsModal.find('.modal-content').html doc
                    # Remove any straggler overlays from rapid re-opening of modal
                    $(".lean-overlay").remove()
                    self.peptideDetailsModal.openModal()

    return TandemGlycopeptideDatabaseSearchResultsController
)()
