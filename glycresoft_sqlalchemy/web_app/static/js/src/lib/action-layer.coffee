class ActionLayerManager extends EventEmitter
    @HOME_LAYER = "home-layer"
    constructor: (container, options) ->
        @container = $(container)
        @options = options
        @layers = {}
        @lastAdded = null
        @layerCounter = 0
        @on 'layer-change', (event) ->
            materialRefresh()

    incLayerCounter: ->
        @layerCounter++
        return @layerCounter

    add: (layer) ->
        @layers[layer.id] = layer
        @container.append layer.container
        @emit "layer-added", "layer": layer
        @lastAdded = layer.id
        return @

    get: (id) ->
        @layers[id]

    getShowingLayer: ->
        result = undefined
        _.forIn @layers, (layer, id) ->
            if layer.showing
                result = layer
            return
        result

    setShowingLayer: (id) ->
        current = @getShowingLayer()
        next = @get(id)
        try
            current.hide()
        try
            next.show()
        catch
            @get(ActionLayerManager.HOME_LAYER).show()
        @emit 'layer-change', 'layer': next

    addLayer: (options, params) ->
        new ActionLayer(this, options, params)
        return @

    removeLayer: (id) ->
        @layers[id].container.remove()
        delete @layers[id]
        return @

    removeCurrentLayer: (next=ActionLayerManager.HOME_LAYER) ->
        current = @getShowingLayer()
        @setShowingLayer next
        current.dispose()
        return @



class ActionLayer
    @actions = {}
    constructor: (manager, options, params) ->
        @manager = manager
        @options = options
        @params = params
        @contentURL = options.contentURL
        if !options.container
            if @params?
                @id = options.name + "-" + manager.incLayerCounter()
            else
                @id = options.name or 'action-layer-' + manager.incLayerCounter()
            @container = $('<div></div>').attr('id', @id)
            @setup()
        else
            @container = $(options.container)
            @id = @container.attr('id')
        @name = options.name or 'layer-' + @id
        @container.attr('data-name', @name).addClass 'container'
        @manager.add this
        @showing = false
        @hide()
        if @manager.getShowingLayer() == undefined
            @show()

    setup: ->
        console.log("Setting up", @)
        if @options.contentURLTemplate?
            @contentURL = @options.contentURLTemplate.format @params
        self = this
        $.get(@contentURL).success (doc) ->
            self.container.hide()
            self.container.html doc
            self.container.find('script').each (i, tag) ->
                tag = $(tag)
                srcURL = tag.attr('src')
                console.log("Setting up script", tag)
                if srcURL != undefined
                    $.getScript srcURL

    show: ->
        @container.fadeIn 100
        @showing = true
        return

    hide: ->
        @container.fadeOut 100
        @showing = false
        return

    dispose: ->
        @container.remove()
        delete @manager.layers[@id]

ActionBook =
    home:
        container: '#home-layer'
        name: 'home-layer'
    addSample:
        contentURL: '/add_sample'
        name: 'add-sample'
    matchSamples:
        contentURL: '/match_samples'
        name: 'match-samples'
    naiveGlycopeptideSearchSpace:
        contentURL: "/glycopeptide_search_space"
        name: "glycopeptide-search-space"
    naiveGlycanSearchSpace:
        contentURL: "/glycan_search_space"
        name: "glycan-search-space"
    viewDatabaseSearchResults:
        contentURLTemplate: "/view_database_search_results/{}"
        name: "view-database-search-results"
