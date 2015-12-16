class ColorManager
    constructor: (obj) ->
        @store = obj
    get: (name) ->
        try
            rgb = @store[name]
            [r, g, b] = rgb
            string = "rgba(#{r * 255},#{g * 255},#{b * 255},0.5)"
            return string
        catch
            @update
    update: (callback) ->
        if !callback?
            callback = () -> {}
        $.getJSON("/api/colors").success (data) =>
            @store = data
            callback(@)



