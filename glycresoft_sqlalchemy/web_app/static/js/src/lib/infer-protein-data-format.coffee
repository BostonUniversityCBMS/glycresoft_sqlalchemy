identifyProteomicsFormat = (file, callback) ->
    isMzidentML = (lines) ->
        for line in lines
            if /mzIdentML/.test line
                return true
        return false

    reader = new FileReader()
    reader.onload = ->
        lines = @result.split("\n")
        console.log lines
        proteomicsFileType = "fasta"
        if isMzidentML(lines)
            proteomicsFileType = "mzIdentML"
        callback(proteomicsFileType)
    reader.readAsText(file.slice(0, 100))


getProteinName = (line) ->
    id = /id="([^"]+)"/.exec(line)
    id = id[1]
    return id.split("_").slice(1).join("_")

getProteinNamesFromMzIdentML = (file, callback, nameCallback) ->
    fr = new FileReader()

    if !nameCallback?
        nameCallback = (name) ->
            console.log(name)
    chunksize = 1024 * 8
    offset = 0
    proteins = {}
    fr.onload = ->
        lines = @result.split("\n")
        for line in lines
            if /<ProteinDetectionHypothesis/i.test line
                name = getProteinName(line)
                if !proteins[name]
                    proteins[name] = true
                    nameCallback(name)

        seek()
    fr.onerror = (error) ->
        console.log(error)
    seek = ->
        if offset >= file.size
            callback Object.keys(proteins)
        else
            fr.readAsText(file.slice(offset, offset + chunksize))
            offset += chunksize / 2
    seek()


class MzIdentMLProteinSelector
    constructor: (file, listContainer) ->
        @fileObject = file
        @container = $(listContainer)
        @initializeContainer()

    initializeContainer: ->
        template = """
        <div class='display-control'>
            <a class='toggle-visible-btn right' data-open="open" style='cursor:hand;'>
                <i class="material-icons">keyboard_arrow_up</i>
            </a>
        </div>
        <div class='hideable-container'>
            <div class='row'>
                <div class='col s4'>
                    <div class='input-field'>
                        <input value='' name="protein-regex" type="text" class="validate protein-regex">
                        <label class="active" for="protein-regex">Protein Pattern</label>
                    </div>
                </div>
                <div class='col s2'>
                    <input type='checkbox' id='select-all-proteins-checkbox' name='select-all-proteins-checkbox'/>
                    <label for='select-all-proteins-checkbox'>Select All</label>
                </div>
            </div>
            <div class='row'>
                <div class='col s8 protein-name-list'>

                </div>
            </div>
        </div>
        """
        @container.html template 
        @hideableContainer = @container.find ".hideable-container"
        @regex = @container.find ".protein-regex"
        @list = @container.find ".protein-name-list"
        @toggleVisible = @container.find ".toggle-visible-btn"
        @selectAllChecker = @container.find '#select-all-proteins-checkbox'
        self = @

        @toggleVisible.click ->
            handle = $ @
            if handle.attr("data-open") == "open"
                self.hideableContainer.hide()
                handle.attr "data-open", "closed"
                handle.html '<i class="material-icons">keyboard_arrow_down</i>'
            else if handle.attr("data-open") == "closed"
                self.hideableContainer.show()
                handle.attr "data-open", "open"
                handle.html '<i class="material-icons">keyboard_arrow_up</i>'

        @regex.change (e) ->
            e.preventDefault()
            pattern = $(@).val()
            self.updateVisibleProteins pattern

        @regex.keydown (e) =>
            if(e.keyCode == 13)
                e.preventDefault();
                @regex.change()
                return false;
        @selectAllChecker.click (e) =>
            if @selectAllChecker.prop("checked")
                @container.find("input[type='checkbox']:visible").prop("checked", true)
            else
                @container.find("input[type='checkbox']:visible").prop("checked", false)
        @load()

    createAddProteinNameToListCallback: ->
        callback = (name) =>
            entryContainer = $("<p></p>").css({"padding-left": 20, "display": 'inline-block', "width": 240}).addClass('input-field protein-name')
            checker = $("<input />").attr("type", "checkbox").attr("name", name).addClass("protein-name-check")
            entryContainer.append checker
            entryContainer.append $("<label></label>").html(name).attr("for", name).click(( -> checker.click()))
            @list.append entryContainer
        return callback

    updateVisibleProteins: (pattern) ->
        regex = new RegExp(pattern, 'i')
        $('.protein-name').each ->
            handle = $ @
            name = handle.find("input").attr("name")
            if regex.test name
                handle.show()
            else
                handle.hide()

    load: ->
        callback = @createAddProteinNameToListCallback()
        getProteinNamesFromMzIdentML(@fileObject, (->), callback)

    getChosenProteins: ->
        return ($(a).attr("name") for a in @container.find("input.protein-name-check:checked"))

    getAllProteins: ->
        return ($(a).attr("name") for a in @container.find("input.protein-name-check"))
