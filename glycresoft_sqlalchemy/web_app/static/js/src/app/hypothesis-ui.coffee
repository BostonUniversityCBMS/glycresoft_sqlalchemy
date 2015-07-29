Application::renderHypothesisListAt = (container)->
    chunks = []
    template = 
    for hypothesis in _.sortBy(_.values(@hypotheses), (o) -> o.name)
        row = $("
    <div data-id=#{hypothesis.id} class=''>
        <span class='handle'>#{hypothesis.name.replace('_', ' ')}</span> <small class='right'>#{if hypothesis.hypothesis_type? then hypothesis.hypothesis_type else '-' }
            <a class='remove-hypothesis mdi-content-clear'></a></small></div>
    ")
        chunks.push row
        row.find(".remove-hypothesis").click (event) -> 
            handle = $ @

    $(container).html chunks

Application.initializers.push ->
    @on "render-hypotheses", =>
        @renderHypothesisListAt ".hypothesis-list"
