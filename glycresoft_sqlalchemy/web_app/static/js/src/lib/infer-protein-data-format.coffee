identifyProteomicsFormat = (file, callback) ->
    isMzidentML = (lines) ->
        for line in lines
            if /mzIdentML/.test line
                return true
        return false

    reader = new FileReader()
    reader.onload = ->
        lines = @result.split("\n")
        console.log lines
        proteomicsFileType = "fasta"
        if isMzidentML(lines)
            proteomicsFileType = "mzIdentML"
        callback(proteomicsFileType)
    reader.readAsText(file.slice(0, 100))


getProteinName = (line) ->
    line.split("_", 2)[1]

getProteinNamesFromMzIdentML = (file, callback) ->
    fr = new FileReader()
    chunksize = 1024 * 8
    offset = 0
    proteins = {}
    fr.onload = ->
        lines = @result.split("\n")
        for line in lines
            if /<ProteinDetectionHypothesis/i.test line
                name = getProteinName(line)
                console.log(name)
                proteins[name] = true
        seek()
    fr.onerror = (error) ->
        console.log(error)
    seek = ->
        if offset >= file.size
            callback Object.keys(proteins)
        else
            fr.readAsText(file.slice(offset, offset + chunksize))
            offset += chunksize / 2
    seek()
