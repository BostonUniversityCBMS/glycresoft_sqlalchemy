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
        @layerStack = []

    incLayerCounter: ->
        @layerCounter++
        return @layerCounter

    add: (layer) ->
        if !layer.options.closeable?
            layer.options.closeable = true
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
        console.log id
        current = @getShowingLayer()
        next = @get(id)
        try
            current.hide()
        try
            next.show()
            i = @findLayer next
            if i != -1
                @layerStack.pop i
            @layerStack.push next
        catch
            @get(ActionLayerManager.HOME_LAYER).show()
        @emit 'layer-change', 'layer': next

    addLayer: (options, params) ->
        if !options.closeable?
            options.closeable = true
        layer = new ActionLayer(this, options, params)
        console.log layer
        if @layerStack.length == 0
            @layerStack.push layer
        return @

    removeLayer: (id) ->
        @layers[id].container.remove()
        i = @findLayer(@layers[id])
        if i != -1
            @layerStack.pop i
        delete @layers[id]
        return @

    removeCurrentLayer: (next=null) ->
        current = @getShowingLayer()
        @layerStack.pop()
        if !next?
            next = @layerStack[@layerStack.length - 1]
        @setShowingLayer next
        current.dispose()
        return @

    findLayer: (targetLayer) ->
        index = -1
        for layer in @layerStack
            index += 1
            if layer.id == targetLayer.id
                return index
        return -1

    setActiveLayerController: (controller) ->
        layer = @getShowingLayer()
        layer.setController(controller)

    getActiveLayerController: ->
        layer = @getShowingLayer()
        return layer.controller


class ActionLayer
    @actions = {}
    constructor: (manager, options, params, method='get') ->
        @manager = manager
        @options = options
        @params = params
        @contentURL = options.contentURL
        @method = method
        if @options.method?
            @method = @options.method
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
        @controller = null

    setup: () ->
        console.log("Setting up", @)
        if @options.contentURLTemplate?
            @contentURL = @options.contentURLTemplate.format @params
        
        callback = (doc) =>
            if not @showing
                @container.hide()
            @options.document = doc
            @container.html doc
            @container.find('script').each (i, tag) ->
                tag = $(tag)
                srcURL = tag.attr('src')
                console.log("Setting up script", tag)
                if srcURL != undefined
                    $.getScript srcURL
            materialRefresh()
            console.log("This layer can be closed? #{@options.closeable}")
            if @options.closeable
                @container.prepend("""
        <div>
            <a class='dismiss-layer mdi-content-clear' onclick='GlycReSoft.removeCurrentLayer()'></a>
        </div>""")
        if @method == "get"
            $.get(@contentURL).success callback
        else if @method == "post"
            $.ajax(@contentURL,
                contentType: "application/json"
                data: JSON.stringify {params: @params, context: @manager.context, settings: @manager.settings}
                success: callback
                type: "POST"
                )
    
    reload: ->
        @container.html @options.document
        @container.find('script').each (i, tag) ->
            tag = $(tag)
            srcURL = tag.attr('src')
            console.log("Setting up script", tag)
            if srcURL != undefined
                $.getScript srcURL
        materialRefresh()
        if @options.closeable
            @container.prepend("""
        <div>
        <a class='dismiss-layer mdi-content-clear' onclick='GlycReSoft.removeCurrentLayer()'></a>
        </div>""")

    setController: (controller) ->
        @controller = controller

    show: ->
        @container.fadeIn 100
        @showing = true

    hide: ->
        @container.fadeOut 100
        @showing = false

    dispose: ->
        @container.remove()
        delete @manager.layers[@id]
