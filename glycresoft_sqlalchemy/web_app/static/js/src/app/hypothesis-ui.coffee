Application::renderHypothesisListAt = (container)->
    chunks = []
    template = ''
    self = @
    for hypothesis in _.sortBy(_.values(@hypotheses), (o) -> o.name)
        if hypothesis.is_decoy
            continue
        row = $("
    <div data-id=#{hypothesis.id} class='list-item clearfix'>
        <span class='handle'>#{hypothesis.name.replace('_', ' ')}</span>
        <small class='right' style='display:inherit'>
            #{if hypothesis.hypothesis_type? then hypothesis.hypothesis_type else '-' }
            <a class='remove-hypothesis mdi-content-clear'></a>
        </small>
    </div>
    ")
        chunks.push row

        row.click (event) ->
            handle = $ @
            hypothesisId = handle.attr("data-id")
            self.addLayer ActionBook.viewHypothesis, {"hypothesis_id": hypothesisId}
            layer = self.lastAdded
            self.setShowingLayer layer
        row.find(".remove-hypothesis").click (event) -> 
            handle = $ @

    $(container).html chunks

Application.initializers.push ->
    @on "render-hypotheses", =>
        @renderHypothesisListAt ".hypothesis-list"
