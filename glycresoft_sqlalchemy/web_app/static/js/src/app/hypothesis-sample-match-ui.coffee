Application::renderHypothesisSampleMatchListAt = (container)->
    chunks = []
    template = 
    for hsm in _.sortBy(_.values(@hypothesisSampleMatches), (o) -> o.id)
        hsm.name = if hsm.name? then hsm.name else "HypothesisSampleMatch:#{hsm.target_hypothesis.name}@#{hsm.sample_run_name}"
        row = $("
    <div data-id=#{hsm.id} class='list-item clearfix'>
        <span class='handle'>#{hsm.id}. #{hsm.name.replace('_', ' ')}</span>
        <small class='right' style='display:inherit'>
            #{hsm.hypothesis_sample_match_type.replace('HypothesisSampleMatch', '')}
            <a class='remove-hsm mdi-content-clear'></a>
        </small>
    </div>
    ")
        chunks.push row
        self = @
        row.click (event) ->
            handle = $(@)
            id = handle.attr('data-id')
            self.addLayer ActionBook.viewDatabaseSearchResults, {hypothesis_sample_match_id: id}
            console.log self.layers
            console.log self.lastAdded
            self.context["hypothesis_sample_match_id"] = id
            self.setShowingLayer self.lastAdded

        row.find(".remove-hsm").click (event) -> 
            handle = $ @
            console.log "Removal of HypothesisSampleMatch is not implemented."

    $(container).html chunks

Application.initializers.push ->
    @on "render-hypothesis-sample-matches", =>
        @renderHypothesisSampleMatchListAt ".hypothesis-sample-match-list"
