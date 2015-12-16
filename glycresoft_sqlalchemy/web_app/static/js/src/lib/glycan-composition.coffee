class GlycanComposition
    constructor: (string) ->
        @__order = []
        @map = {}
        @parse(string)

    parse: (string) ->
        parts = string.slice(1, -1).split("; ")
        for part in parts
            [name, number] = part.split(":")
            @__order.push(name)
            @map[name] = parseInt(number)

    format: (colorSource) ->
        parts = []
        for name, number of @map
            if name == '__order'
                continue
            color = colorSource.get(name)
            template = """<span class='monosaccharide-name' style='background-color:#{color}; padding: 2px;border-radius:2px;'>#{name} #{number}</span>"""
            parts.push template
        return parts.join ' '
