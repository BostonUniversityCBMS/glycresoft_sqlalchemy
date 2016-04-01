module.exports = (grunt) -> 
    grunt.initConfig
        pkg: grunt.file.readJSON("package.json")
        copy:
            js:
                flatten: true
                expand: true
                src: 'static/js/dist/*'
                dest: "static/dist/js/"
            css:
                flatten: true
                expand: true
                src: "static/css/dist/*"
                dest: "static/dist/css/"
            icons:
                flatten: true
                expand: true
                src: "static/font/material-design-icons/*"
                dest: "static/dist/font/material-design-icons/"
            icons2:
                flatten: true
                expand: true
                src: "static/font/material-design-icons/*"
                dest: "static/dist/css/material-design-icons/"
            font:
                flatten: true
                expand: true
                src: "static/font/roboto/*"
                dest: "static/dist/font/roboto/"
        concat:
            options:
                separator: "\n"
            app:
                src: ["static/js/build/app/**.js"]
                dest: "static/js/dist/app-bundle.js"
            lib:
                src: ["static/js/build/lib/**.js"]
                dest: "static/js/dist/lib-bundle.js"
            vendor:
                src: ["static/js/vendor/**.js"]
                dest: "static/js/dist/vendor-bundle.js"
            css:
                src: ["static/css/*.css"]
                dest: "static/css/dist/bundle.css"
                separator: "\n"
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
    grunt.loadNpmTasks "grunt-contrib-copy"


    grunt.registerTask "default", ["coffee", "concat", "copy"]

    return grunt
