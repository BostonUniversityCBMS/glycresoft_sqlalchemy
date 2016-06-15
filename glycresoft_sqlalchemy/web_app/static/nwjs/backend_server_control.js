var GetNextPort = process.mainModule.exports.GetNextPort
var TerminateServer = process.mainModule.exports.TerminateServer

BackendServerControl = (function(){
var http = require('http')
var gui = require('nw.gui')
var child_process = require("child_process")

WINDOW_OPTIONS = {
    "title": "GlycReSoft",
    "toolbar": true,
    "icon": "static/nwjs/logo.png"
}


function _getServerExecutable(){
    var argIx = gui.App.argv.indexOf("--executable")
    if(argIx == -1){
        return "glycresoft-report web"
    } else {
        return gui.App.argv[argIx + 1]
    }
}


var EXECUTABLE = _getServerExecutable()


function BackendServerControl(project, options){
    options = options === undefined ? {} : options;
    this.project = project
    this.port = options.port === undefined ? GetNextPort() : options.port
    this.host = options.host === undefined ? "localhost" : options.host
    this.protocol = options.protocol === undefined ? "http:" : options.protocol
    this.terminateCallback = options.callback === undefined ? function(){} : options.callback
    this.url = "http://localhost:" + this.port
    this.process = null
    this.window = null
}

BackendServerControl.EXECUTABLE = EXECUTABLE
BackendServerControl.prototype.EXECUTABLE = EXECUTABLE

BackendServerControl.prototype.constructServerProcessCall = function(){
    return this.EXECUTABLE + " " + this.project.getStorePath() + " --port " + this.port
}

BackendServerControl.prototype.launchServer = function(){
    child = child_process.exec(this.constructServerProcessCall())
    child.stdout.on("data", function(){
        console.log("stdout", arguments)
        console.log(arguments[0])
    })
    child.stderr.on("data", function(){
        console.log("stderr", arguments)
        console.log(arguments[0])
    })
    this.process = child
}


BackendServerControl.prototype.configureTerminationBehavior = function(){
    var self = this
    self.process.on("exit", function(){
        console.log("Exited!", arguments)
        self.terminateCallback(self)
        self.window.close(true)
    })
    self.window.on("close", function(){
        try {
            self.terminateServer()
            self.process.kill()
        } catch(error){
            console.log("Could not terminate server", error)
        }
        self.terminateCallback(self)
        self.window.close(true)
    })

}

BackendServerControl.prototype.navigateOnReady = function(count, callback){
    var url = this.url
    var self = this
    count = count === undefined ? 1 : count + 1;
    if(count > 600){
        throw new Error("Server Not Ready After " + count + " Tries")
    }
    console.log("Calling navigateOnReady with", self, count)
    http.get(self.url, function(response){
        var retry = false
        if(response.statusCode == 200){
            console.log(self.url)
            try{
                self.window = gui.Window.open(self.url, WINDOW_OPTIONS)
                self.window.maximize()
                if(callback !== undefined){
                    callback()
                }                
            } catch(error){
                retry = true
                console.log(error)
            }
        } else {
            retry = true
        }
        if(retry){
            self.navigateOnReady(count, callback)
        }
    }).on('error', function(e) {
      setTimeout(function(){self.navigateOnReady(count, callback)}, 150)
    });
}

BackendServerControl.prototype.terminateServer = function(){
    http.request({host:this.host, "port": this.port, protocol: this.protocol, path: "/internal/shutdown", method: "POST"}).end()
}

BackendServerControl.launch = function(project, options){
    var server = new BackendServerControl(project, options)
    server.launchServer()
    console.log("Server Launched")
    server.navigateOnReady(0, function(){server.configureTerminationBehavior()})
    return server
}

return BackendServerControl
})()
