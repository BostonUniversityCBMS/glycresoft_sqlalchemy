class Application extends ActionLayerManager
    constructor: (options={}) ->
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

        self = this
        @eventStream = new EventSource('/stream')
        
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
        @handleMessage "new-sample", (data) =>
            console.log(data)       

    runInitializers: ->
        console.log Application, Application.initializers
        for initializer in Application.initializers
            console.log initializer
            initializer.apply this, null 

    updateSettings: ->
        $.post('/internal/update_settings', @settings).success((data) ->
            @settings = data
        ).error (err) ->
            console.log err

    updateTaskList: ->
        taskListContainer = @sideNav.find('.task-list-container ul')

        clickTask = (event) ->
            handle = $(this)
            state = handle.attr('data-status')
            id = handle.attr('data-id')
            if state == 'finished'
                console.log self.tasks[id]
                delete self.tasks[id]
                handle.fadeOut()
                handle.remove()
            return

        taskListContainer.html _.map(@tasks, renderTask).join('')
        self = this
        taskListContainer.find('li').click clickTask

    handleMessage: (messageType, handler) ->
        @eventStream.addEventListener messageType, (event) ->
            data = JSON.parse(event.data)
            handler(data)

    @initializers = [
        ->
            console.log this
        ->
            self = this
            $ ->
                self.container = $(self.options.actionContainer)
                console.log self.options.actionContainer
                self.sideNav = $('.side-nav')
                self.addLayer ActionBook.home
    ]

renderTask = (task) ->
    '<li data-id=\'{id}\' data-status=\'{status}\'><b>{name}</b> ({status})</li>'.format task


window.GlycReSoft = new Application(options={actionContainer: ".action-layer-container"})

$(() ->
    console.log("updating Application")
    GlycReSoft.runInitializers()
    GlycReSoft.updateSettings()
    GlycReSoft.updateTaskList())
