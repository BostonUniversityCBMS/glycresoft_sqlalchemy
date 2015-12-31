
class MonosaccharideFilter
    constructor: (parent, residueNames, rules) ->
        if !rules?
            if !GlycReSoft.settings.monosaccharide_filters?
                GlycReSoft.settings.monosaccharide_filters = {}
            rules = GlycReSoft.settings.monosaccharide_filters
        @container = $("<div></div>").addClass("row")
        $(parent).append(@container)
        @residueNames = residueNames
        @rules = rules

    makeFilterWidget: (residue) ->
        rule = @rules[residue]
        if !rule?
            rule = {
                minimum: 0
                maximum: 10
                include: true
            }
            @rules[residue] = rule
        residue.name = residue
        residue.sanitizeName = sanitizeName = residue.replace(/[\(\),]/g, "_")
        template = """
            <span class="col s2 monosaccharide-filter" data-name='#{residue}'>
                <p style='margin: 0px; margin-bottom: -10px;'>
                    <input type="checkbox" id="#{sanitizeName}_include" name="#{sanitizeName}_include"/>
                    <label for="#{sanitizeName}_include"><b>#{residue}</b></label>
                </p>
                <p style='margin-top: 0px; margin-bottom: 0px;'>
                    <input id="#{sanitizeName}_min" type="number" placeholder="Minimum #{residue}" style='width: 45px;' min="0"
                           value="#{rule.minimum}" max="#{rule.maximum}" name="#{sanitizeName}_min"/> : 
                    <input id="#{sanitizeName}_max" type="number" placeholder="Maximum #{residue}" style='width: 45px;' min="0"
                           value="#{rule.maximum}" name="#{sanitizeName}_max"/>
                </p>
            </span>
            """
        self = @
        rendered = $(template)
        rendered.find("##{sanitizeName}_min").change ->
            rule.minimum = parseInt($(@).val())
            self.changed()
        rendered.find("##{sanitizeName}_max").change ->
            rule.maximum = parseInt($(@).val())
            self.changed()
        rendered.find("##{sanitizeName}_include").prop("checked", rule.include).click ->
            rule.include = $(@).prop("checked")
            self.changed()

        return rendered

    render: ->
        for residue in @residueNames
            widget = @makeFilterWidget(residue)
            @container.append(widget)

    changed: _.debounce (-> GlycReSoft.emit("update_settings")), 1000 
