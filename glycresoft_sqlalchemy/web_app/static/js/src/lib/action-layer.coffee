class ActionLayerManager extends EventEmitter
    @HOME_LAYER = "home-layer"
    constructor: (container, options) ->
        @container = $(container)
        @options = options
        @layers = {}
        @layerCounter = 0
        @on 'layer-change', (event) ->
            materialRefresh()
            return
        return

    incLayerCounter: ->
        @layerCounter++
        @layerCounter

    add: (layer) ->
        @layers[layer.id] = layer
        @container.append layer.container
        return

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

    addLayer: (options) ->
        new ActionLayer(this, options)

    removeLayer: (id) ->
        @layers[id].container.remove()
        delete @layers[id]

class ActionLayer
    constructor: (manager, options) ->
        @manager = manager
        @options = options
        @contentURL = options.contentURL
        if !options.container
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
        self = this
        $.get(@contentURL).success (doc) ->
            self.container.hide()
            self.container.html doc
            self.container.find('script').each (i, tag) ->
                tag = $(tag)
                srcURL = tag.attr('src')
                if srcURL != undefined
                    $.getScript srcURL
                else
                    eval tag.text()

    show: ->
        @container.fadeIn 100
        @showing = true
        return

    hide: ->
        @container.fadeOut 100
        @showing = false
        return
