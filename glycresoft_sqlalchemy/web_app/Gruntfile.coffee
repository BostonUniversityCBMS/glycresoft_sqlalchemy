module.exports = (grunt) -> 
    grunt.initConfig
        pkg: grunt.file.readJSON("package.json")
        concat:
            options:
                separator: ";"
            app:
                src: ["static/js/build/app/**.js"]
                dest: "static/js/dist/app-bundle.js"
            lib:
                src: ["static/js/build/lib/**.js"]
                dest: "static/js/dist/lib-bundle.js"
            vendor:
                src: ["static/js/vendor/**.js"]
                dest: "static/js/dist/vendor-bundle.js"
        coffee:
            app:
                options:
                    bare: true
                    sourceMap: true
                expand: true,
                cwd: "static/js/src/app/"
                src: ["*.coffee"]
                dest: "static/js/build/app/"
                ext: ".js"
            lib:
                options:
                    bare: true
                    sourceMap: true
                expand: true,
                cwd: "static/js/src/lib/"
                src: ["*.coffee"]
                dest: "static/js/build/lib/"
                ext: ".js"

    grunt.loadNpmTasks 'grunt-contrib-concat'
    grunt.loadNpmTasks 'grunt-contrib-coffee'


    grunt.registerTask "default", ["coffee", "concat"]

    return grunt
