(function(){
    var child_process = require("child_process")
    var http = require('http')

    var basePort = 5000
    var portOffset = 0;
    function GetNextPort(){
        var port = basePort + portOffset;
        portOffset += 1;
        return port
    }

    function TerminateServer(port){
        http.request({host:"localhost", "port": port, path: "/internal/shutdown", method: "POST"}).end()
    }

    exports.GetNextPort = GetNextPort
    exports.TerminateServer = TerminateServer
})();
