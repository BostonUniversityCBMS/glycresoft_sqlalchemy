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
    constructor: (manager, options, params, method='get') ->
        @manager = manager
        @options = options
        @params = params
        @contentURL = options.contentURL
        @method = method
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
        
        callback = (doc) =>
            if not @showing
                @container.hide()
            @container.html doc
            @container.find('script').each (i, tag) ->
                tag = $(tag)
                srcURL = tag.attr('src')
                console.log("Setting up script", tag)
                if srcURL != undefined
                    $.getScript srcURL
            materialRefresh()
            @container.prepend("""
    <div>
        <a class='dismiss-layer mdi-content-clear' onclick='GlycReSoft.removeCurrentLayer()'></a>
    </div>""")
        if @method == "get"
            $.get(@contentURL).success callback
        else if @method == "post"
            $.post(@contentURL, @params).success callback
        

    show: ->
        @container.fadeIn 100
        @showing = true

    hide: ->
        @container.fadeOut 100
        @showing = false

    dispose: ->
        @container.remove()
        delete @manager.layers[@id]
