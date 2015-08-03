Application::renderSampleListAt = (container)->
    chunks = []
    template = 
    for sample in _.sortBy(_.values(@samples), (o) -> o.name)
        row = $("
    <div data-name=#{sample.name}>
        <span class='handle'>#{sample.name.replace('_', ' ')}</span> <small class='right'>#{sample.sample_type}
            <a class='remove-sample mdi-content-clear'></a></small></div>
    ")
        chunks.push row
        row.find(".remove-sample").click (event) -> 
            handle = $ @
            console.log handle

    $(container).html chunks

Application.initializers.push ->
    @on "render-samples", =>
        @renderSampleListAt ".sample-list"
