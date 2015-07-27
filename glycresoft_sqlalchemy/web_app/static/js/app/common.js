var Application = function(){
    this.version = [0, 0, 1]
    this.context = {}
    this.settings = {}
    this.tasks = {}
    this.sideNav = $(".side-nav")

    var self = this
    this.eventStream = new EventSource('/stream');
    this.eventStream.addEventListener("update", function(event){
        Materialize.toast(event.data.replace(/"/g, ''), 4000)
    })
    this.eventStream.addEventListener("task-queued", function(event){
        data = JSON.parse(event.data)
        self.tasks[data.id] = {'id': data.id, "name": data.name, "status": "queued"}
        self.updateTaskList()
    })    
    this.eventStream.addEventListener("task-start", function(event){
        data = JSON.parse(event.data)
        self.tasks[data.id] = {'id': data.id, "name": data.name, "status": "running"}
        self.updateTaskList()
    })
    this.eventStream.addEventListener("task-complete", function(event){
        data = JSON.parse(event.data)
        try{
            self.tasks[data.id].status = "finished"                    
        } catch(err){
            self.tasks[data.id] = {'id': data.id, "name": data.name, "status": "finished"}
        }
        self.updateTaskList()
    })

    for(var i = 0; i < Application.initializers.length; i++){
        var initializer = Application.initializers[i]
        initializer.apply(this, null)
    }
}

_.extend(Application.prototype, EventEmitter.prototype)


Application.prototype.updateSettings = function(){
    $.post('/internal/update_settings', this.settings).success(function(data){
        this.settings = data
    }).error(function(err){
        console.log(err)
    })
}


function renderTask(task){
    return "<li data-id='{id}' data-status='{status}'><b>{name}</b> ({status})</li>".format(task)
}


Application.prototype.updateTaskList = function(){
    var taskListContainer = this.sideNav.find(".task-list-container ul")
    taskListContainer.html(_.map(this.tasks, renderTask).join(''))
    var self = this

    function clickTask(event){
        var handle = $(this)
        var state = handle.attr('data-status')    
        var id = handle.attr('data-id')
        if (state === 'finished'){
            console.log(self.tasks[id])
            delete self.tasks[id]
            handle.fadeOut()
            handle.remove()
        }
    }

    taskListContainer.find("li").click(clickTask)
}

Application.initializers = [
    function(){console.log(this)},
    function(){
        var self = this
        $(function(){
            self.sideNav = $(".side-nav")
        })
    }
    ]

window.GlycReSoft = new Application()

$(function(){
    GlycReSoft.updateSettings()
    GlycReSoft.updateTaskList()
})
