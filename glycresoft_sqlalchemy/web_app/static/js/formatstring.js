
(function(){
    /*
    Implements {named} replacements similar to the simple format() method of strings from Python
     */
    String.prototype.format = function(){
        var data = arguments
        var i = 0;
        var keys = Object.keys(arguments);
        if(arguments.length === 1 && (typeof arguments[0] === 'object')){
            data = arguments[0]
            keys = Object.keys(arguments);
        }
        res = this.replace(/\{([^\}]*)\}/g, function(placeholder, name, position){
            if (name === ""){
                name = keys[i];
                i++;
            }
            try{
                return JSON.stringify(data[name]).slice(1, -1)
            } catch(err){
                console.log(err, name, data)
                return undefined
            }
        })
        return res
    }
})()
