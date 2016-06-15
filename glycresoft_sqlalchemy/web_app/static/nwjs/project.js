parts = (function(){

var gui = require("nw.gui")
var rimraf = require("rimraf")
var fs = require('fs')

var WINDOW = gui.Window.get()
var PROJECTS_KEY = "projects"
var PROJECT_FILE = "store.glycresoftdb"
var VERSION = 0.1


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
            if(project_i.path === project.path){
                console.log(project_i)
                continue;
            }
            filter.push(project_i)
        }
        localforage.setItem(PROJECTS_KEY, filter)
    })
}


function makeNewProjectDirectory(path, name){
    name = name.replace(/\s/g, '_')
    path = [path, name].join('/')
    fs.mkdirSync(path)
    return path
}


function MakeProjectFromDOM(){
    var path = $("input[name='project-location']").val()
    var name = $("#project_name").val()
    path = makeNewProjectDirectory(path, name)
    return new Project(name, path)
}

function ProjectSelectionWindow(){
    this.controllers = []
    this.projects = []
    this.window = gui.Window.get()
    this.updateProjectDisplay()
    var self = this
    $("#create-project-btn").click(function(){self.createProject()})
    $("#delete-existing-btn").click(function(){self.deleteProject()})
    $("#load-existing-btn").click(function(){self.openProject()})

    process.on("exit", function(){
        for(var i=0;i<self.controllers.length;i++){
            var controller = self.controllers[i]
            controller.terminateServer()
        }
    })

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
    var self = this
    AddProjectToLocalStorage(project).then(function(){self.updateProjectDisplay();})
    this.controllers.push(BackendServerControl.launch(
        project, {callback: this.dropBackendController.bind(this)}))
}

ProjectSelectionWindow.prototype.deleteProject = function(){
    var selectProjectTag = $("select#existing-project"); 
    var project = this.projects[selectProjectTag.val()]
    var self = this
    console.log(project)
    _RemoveProject(project)
    rimraf(project.path, function(){
        self.updateProjectDisplay()
    })
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
