parts = (function(){

var gui = require("nw.gui")

var WINDOW = gui.Window.get()
var PROJECTS_KEY = "projects"
var PROJECT_FILE = "store.glycresoftdb"
var VERSION = 0.1

var PROJECT_SCHEMA = {
    name: String,
    path: String,
    version: Number
}


function Project(name, path, version){
    if(name instanceof Object){
        this.name = name.name;
        this.path = name.path;
        this.version = name.version
    }
    else {
        this.name = name;
        this.path = path;
        this.version = version === undefined ? VERSION : version;
    }
}

Project.fromObject = function(obj){
    return new Project(obj)
}

Project.prototype.getStorePath = function(){
    return this.path + "/" + PROJECT_FILE
}


function AddProjectToLocalStorage(project){
    return localforage.getItem(PROJECTS_KEY).then(function(value){
        console.log(arguments)
        if((value === undefined) || (value === null)){
            value = [];
        }
        value.push(project)
        return localforage.setItem(PROJECTS_KEY, value)
    })
}

function LoadAllProjects(){
    return localforage.getItem(PROJECTS_KEY).then(function(values){
        return Promise.resolve(values.map(Project.fromObject))
    })
}


function _RemoveAllProjects(){
    return localforage.setItem(PROJECTS_KEY, []);
}


function _RemoveProject(project){
    console.log(arguments)
    return localforage.getItem(PROJECTS_KEY).then(function(value){
        if(value === undefined || value === null){
            value = [];
        }
        filter = []
        for(var i = 0; i < value.length; i++){
            var project_i = value[i];
            if(project_i.location === project.location){
                continue;
            }
            filter.push(project_i)
        }
        localforage.setItem(PROJECTS_KEY, filter)
    })
}

function MakeProjectFromDOM(){
    var path = $("input[name='project-location']").val()
    var name = $("#project_name").val()
    return new Project(name, path)
}

function ProjectSelectionWindow(){
    this.controllers = []
    this.projects = []
    this.window = gui.Window.get()
    this.updateProjectDisplay()
    self = this
    $("#create-project-btn").click(function(){self.createProject()})
    $("#load-existing-btn").click(function(){self.openProject()})

    self.window.on("close", function(){
        if(self.controllers.length == 0){
            self.window.close(true);
        } else {
            self.window.hide()
        }
    })

}

ProjectSelectionWindow.prototype.dropBackendController = function(server){
    console.log("dropBackendController", server)
    var ix = -1;
    for(var i = 0; i < this.controllers.length; i++){
        if(this.controllers[i].port === server.port){
            ix = i
        }
    }
    if(ix != -1){
        this.controllers.pop(ix);
    }
    this.window.close()
}

ProjectSelectionWindow.prototype.openProject = function(){
    var selectProjectTag = $("select#existing-project")
    var project = this.projects[selectProjectTag.val()]
    console.log(project)
    this.controllers.push(BackendServerControl.launch(
        project, {callback: this.dropBackendController.bind(this)}))
}

ProjectSelectionWindow.prototype.updateProjectDisplay = function(){
    var self = this
    LoadAllProjects().then(function(projects){
        var existingContainer = $("#load-existing-project-container");
        if(projects == null || projects.length == 0){
            existingContainer.hide()
            return;
        }
        var selectProjectTag = $("select#existing-project");
        self.projects = projects
        console.log(projects)
        selectProjectTag.empty()
        for(var i = 0; i < projects.length; i++){
            var project = projects[i]
            var displayName = project.name === undefined ? project.path : project.name;
            if(project.path === undefined){
                continue;
            }
            var optionTag = $("<option></option>").text(displayName).attr("value", i)
            selectProjectTag.append(optionTag);
        }
        existingContainer.show();
    })
}

ProjectSelectionWindow.prototype.createProject = function(){
    var project = MakeProjectFromDOM();
    AddProjectToLocalStorage(project);
    this.updateProjectDisplay();
    this.controllers.push(BackendServerControl.launch(
        project, {callback: this.dropBackendController.bind(this)}))
}


global.ProjectSelectionWindow = ProjectSelectionWindow
global.Project = Project
return [Project, ProjectSelectionWindow]
})()

Project = parts[0]
ProjectSelectionWindow = parts[1]

var Controller = null

$(function(){
    Controller = new ProjectSelectionWindow()
})
