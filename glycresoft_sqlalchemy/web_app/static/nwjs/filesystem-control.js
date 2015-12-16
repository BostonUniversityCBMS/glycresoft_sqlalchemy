(function(){
    var fs = require('fs')

    FileTypeTitles = {
        "txt": "Text",
        "fa": "Fasta",
        "fasta": "Fasta",
        "fas": "Fasta",
        "mzid": "mzIdentML",
        "yaml": "YAML",
        "csv": "CSV",
        "db": "Data Store"
    }

    function fileTypeToTitle(fileType){
        var type = FileTypeTitles[fileType]
        if(type === undefined){
            return fileType
        }
        return type + " File"
    }

    function FileElement(filePath){
        this.fullPath = filePath
        this.name = filePath.split("/").slice(-1)[0]
        this.type = this.name.split('.').slice(-1)[0]
    }

    FileElement.prototype.transform = function(container) {
        container.find(".element-icon").append('<i class="material-icons">insert_drive_file</i>')
    }

    function DirectoryElement(dirPath){
        this.fullPath = dirPath
        this.name = dirPath.split("/").slice(-1)[0]
        this.type = "Directory"
    }

    DirectoryElement.prototype.transform = function(container) {
        container.addClass("directory-entry")
        container.find(".element-icon").append('<i class="material-icons">folder</i>')
    }

    function ProjectElement(projectPath){
        this.fullPath = projectPath
        this.name = projectPath.split("/").slice(-1)[0]
        this.type = "GlycReSoft Project"    
    }

    ProjectElement.prototype.transform = function(container) {
        container.addClass("project-entry")
        container.find(".element-icon").append('<i class="material-icons">open_in_browser</i>')
    }

    function detectProject(path){
        var contents = fs.readdirSync(path)
        for(var i = 0; i < contents.length; i++){
            var fileName = contents[i]
            if (fileName == "store.glycresoftdb"){
                return true
            }
        }
        return false
    }

    var PROJECT = "PROJECT"
    var DIRECTORY = "DIRECTORY"
    var FILE = "FILE"

    var TYPE_MAP = {
        PROJECT: ProjectElement,
        DIRECTORY: DirectoryElement,
        FILE: FileElement
    }

    var PATH_CACHE = {}

    function classify(path){
        var stat = fs.statSync(path)
        if(stat.isDirectory()){
            if (detectProject(path)){
                return PROJECT
            } else {
                return DIRECTORY
            }
        } else {
            return FILE
        }
    }

    function constructElementsForPath(path){
        var contents = fs.readdirSync(path)
        var listingResults = []
        for(var i = 0; i < contents.length; i++){
            var file = contents[i];
            var path = fs.realpathSync(file, PATH_CACHE)
            var type = classify(path)
            console.log(type)
            var element = new TYPE_MAP[type](path)
            listingResults.push(element)
        }
        return listingResults
    }

    if (module.exports){    
        module.exports = {
            "constructElementsForPath": constructElementsForPath,
            "detectProject": detectProject,
            "ProjectElement": ProjectElement,
            "DirectoryElement": DirectoryElement,
            "FileElement": FileElement,
            "FileTypeTitles": FileTypeTitles,
        }
    } else {
        var _exported = {
            "constructElementsForPath": constructElementsForPath,
            "detectProject": detectProject,
            "ProjectElement": ProjectElement,
            "DirectoryElement": DirectoryElement,
            "FileElement": FileElement,
            "FileTypeTitles": FileTypeTitles,
        }
        for(var name in _exported){
            window[name] = _exported[name]
        }
    }
})()
