class Application extends ActionLayerManager
    constructor: (options={}) ->
        console.log "Instantiating Application", this

        super options.actionContainer, options

        @version = [
            0
            0
            1
        ]
        @context = {}
        @settings = {}
        @tasks = {}
        @sideNav = $('.side-nav')
        @colors = new ColorManager()
        self = this
        
        @connectEventSource()
        
        @handleMessage 'update', (data) =>
            Materialize.toast data.replace(/"/g, ''), 4000
            return

        @handleMessage 'task-queued', (data) =>
            self.tasks[data.id] =
                'id': data.id
                'name': data.name
                'status': 'queued'
            self.updateTaskList()
            return
        @handleMessage 'task-start', (data) =>
            self.tasks[data.id] =
                'id': data.id
                'name': data.name
                'status': 'running'
            self.updateTaskList()
            return
        @handleMessage 'task-error', (data) =>
            task = self.tasks[data.id]
            task.status = 'error'
            self.updateTaskList()
            return
        @handleMessage 'task-complete', (data) =>
            try
                self.tasks[data.id].status = 'finished'
            catch err
                self.tasks[data.id] =
                    'id': data.id
                    'name': data.name
                    'status': 'finished'
            self.updateTaskList()
            return
        @handleMessage 'new-sample', (data) =>
            @samples[data.id] = data
            @emit "render-samples"
        @handleMessage 'new-hypothesis', (data) =>
            @hypotheses[data.id] = data
            @emit "render-hypotheses"
        @handleMessage 'new-hypothesis-sample-match', (data) =>
            @hypothesisSampleMatches[data.id] = data
            @emit "render-hypothesis-sample-matches"

        @on "layer-change", (data) =>
            @colors.update()

    connectEventSource: ->
        console.log "Establishing EventSource connection"
        @eventStream = new EventSource('/stream')

    runInitializers: ->
        for initializer in Application.initializers
            initializer.apply this, null 

    updateSettings: (payload={}) ->
        $.post('/preferences', payload).success((data) =>
            for k, v of data
                @settings[k] = v
            @emit("update_settings")
        ).error (err) ->
            console.log "error in updateSettings", err, arguments

    updateTaskList: ->
        taskListContainer = @sideNav.find('.task-list-container ul')

        clickTask = (event) ->
            handle = $(this)
            state = handle.attr('data-status')
            id = handle.attr('data-id')
            if state == 'finished'
                delete self.tasks[id]
                handle.fadeOut()
                handle.remove()
            return
        self = @
        doubleClickTask = (event) ->
            handle = $(this)    
            id = handle.attr('data-id')
            $.get "/internal/log/" + id, (message) => self.displayMessageModal(message)

        taskListContainer.html _.map(@tasks, renderTask).join('')
        taskListContainer.find('li').map (i, li) -> contextMenu li, {"View Log": doubleClickTask}
        taskListContainer.find('li').click clickTask
        taskListContainer.find("li").dblclick doubleClickTask

    handleMessage: (messageType, handler) ->
        @eventStream.addEventListener messageType, (event) ->
            data = JSON.parse(event.data)
            handler(data)

    @initializers = [
        ->
            self = @
            $ ->
                self.container = $(self.options.actionContainer)
                self.sideNav = $('.side-nav')
                self.addLayer ActionBook.home
                $("#run-matching").click (event) ->
                    $(".lean-overlay").remove()
                    setupAjaxForm "/ms1_or_ms2_choice?ms1_choice=peakGroupingMatchSamples&ms2_choice=tandemMatchSamples",
                                  "#message-modal"

                $("#build-glycan-search-space").click (event) ->
                    self.addLayer ActionBook.naiveGlycanSearchSpace
                    self.setShowingLayer self.lastAdded
                $("#build-glycopeptide-search-space").click (event) ->
                    self.addLayer ActionBook.naiveGlycopeptideSearchSpace
                    self.setShowingLayer self.lastAdded
        ->
            @loadData()
        ->
            @handleMessage "files-to-download", (data) =>
                for file in data.files
                    @downloadFile(file)
        ->
            @on "update_settings", =>
                layer = @getShowingLayer()
                if layer.name != ActionBook.home.name
                    console.log("Updated Settings, Current Layer:", layer.name)
                    layer.setup()
    ]

    loadData: ->
        DataSource.hypotheses (d) => 
            console.log('hypothesis', d)
            @hypotheses = d
            @emit "render-hypotheses"
        DataSource.samples (d) =>
            console.log('samples', d)
            @samples = d
            @emit "render-samples"
        DataSource.hypothesisSampleMatches (d) =>
            console.log('hypothesisSampleMatches', d)
            @hypothesisSampleMatches = d
            @emit "render-hypothesis-sample-matches"
        DataSource.tasks (d) =>
            console.log('tasks', d)
            @tasks = d
            @updateTaskList()
        @colors.update()

    downloadFile: (filePath) ->
        window.location = "/internal/file_download/" + btoa(filePath)

    displayMessageModal: (message) ->
        container = $("#message-modal")
        container.find('.modal-content').html message
        container.openModal()

    ajaxWithContext: (url, options) ->
        if !options?
            options = {data:{}}
        data = options.data
        data['settings'] = @settings
        data['context'] = @context
        options.method = "POST"
        options.data = JSON.stringify(data)
        options.contentType = "application/json"
        return $.ajax(url, options)



renderTask = (task) ->
    '<li data-id=\'{id}\' data-status=\'{status}\'><b>{name}</b> ({status})</li>'.format task


$(() ->
    window.GlycReSoft = new Application(options={actionContainer: ".action-layer-container"})
    console.log("updating Application")
    GlycReSoft.runInitializers()
    GlycReSoft.updateSettings()
    GlycReSoft.updateTaskList())
