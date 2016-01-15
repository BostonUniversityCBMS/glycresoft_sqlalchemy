_formatModification = (modification, colorSource, long=false) ->
        color = colorSource.get(modification)
        content = escape(if long then modification else modification[0])
        return "(<span class='modification-chip' style='background-color:#{color};padding-left:1px;padding-right:2px;border-radius:2px;'
           title='#{escape(modification)}' data-modification=#{escape(modification)}>#{content}</span>)"


class PeptideSequencePosition
    formatModification = (modification, colorSource, long=false) ->
        color = colorSource.get(modification)
        content = if long then modification else modification[0]
        return "(<span class='modification-chip' style='background-color:#{color};padding-left:1px;padding-right:2px;border-radius:2px;'
           title='#{modification}' data-modification=#{modification}>#{content}</span>)"

    constructor: (residue, modifications=[]) ->
        @residue = residue
        @modifications = modifications

    format: (colorSource) ->
        modifications = @modifications.map((modification) -> formatModification(modification, colorSource)).join ''
        return @residue + modifications

class PeptideSequence
    parser = (sequence) ->
        state = "start"  # [start, nTerm, aa, mod, cTerm]
        nTerm = nTerm or "H"
        cTerm = cTerm or "OH"
        mods = []
        chunks = []
        glycan = ""
        currentAA = ""
        currentMod = ""
        currentMods = []
        parenLevel = 0
        i = 0
        while i < sequence.length
            nextChr = sequence[i]
            # Transition from aa to mod when encountering the start of a modification
            # internal to the sequence
            if nextChr == "("
                if state == "aa"
                    state = "mod"
                    parenLevel += 1
                # Transition to nTerm when starting on an open parenthesis
                else if state == "start"
                    state = "nTerm"
                    parenLevel += 1
                else
                    parenLevel += 1
                    if not (state in ["nTerm", "cTerm"] and parenLevel == 1)
                        currentMod += nextChr
            else if nextChr == ")"
                if state == "aa"
                    throw new Exception(
                        "Invalid Sequence. ) found outside of modification, Position {0}. {1}".format(i, sequence))
                else
                    parenLevel -= 1
                    if parenLevel == 0
                        mods.push(currentMod)
                        currentMods.push(currentMod)
                        if state == "mod"
                            state = 'aa'
                            if currentAA == ""
                                chunks.slice(-1)[1] = chunks.slice(-1)[1].concat currentMods
                            else 
                                chunks.push([currentAA, currentMods])
                        else if state == "nTerm"
                            if sequence[i+1] != "-"
                                throw new Exception("Malformed N-terminus for " + sequence)
                            # Only one modification on termini
                            nTerm = currentMod
                            state = "aa"
                            # Jump ahead past - into the amino acid sequence
                            i += 1
                        else if state == "cTerm"
                            # Only one modification on termini
                            cTerm = currentMod

                        currentMods = []
                        currentMod = ""
                        currentAA = ""
                    else
                        currentMod += nextChr

            else if nextChr == "|"
                if state == "aa"
                    throw new Exception(
                        "Invalid Sequence. | found outside of modification")
                else
                    currentMods.push(currentMod)
                    mods.push(currentMod)
                    currentMod = ""
            else if nextChr == "{"
                if (state == 'aa' or (state == "cTerm" and parenLevel == 0))
                    glycan = sequence.slice(i)
                    break
                else
                    currentMod += nextChr

            else if nextChr == "-"
                if state == "aa"
                    state = "cTerm"
                    if(currentAA != "")
                        currentMods.push(currentMod)
                        chunks.push([currentAA, currentMods])
                        currentMod = ""
                        currentMods = []
                        currentAA = ""
                else
                    currentMod += nextChr
            else if state == "start"
                state = "aa"
                currentAA = nextChr
            else if state == "aa"
                if(currentAA != "")
                    currentMods.push(currentMod)
                    chunks.push([currentAA, currentMods])
                    currentMod = ""
                    currentMods = []
                    currentAA = ""
                currentAA = nextChr
            else if state in ["nTerm", "mod", "cTerm"]
                currentMod += nextChr
            else
                throw new Exception(
                    "Unknown Tokenizer State", currentAA, currentMod, i, nextChr)
            i += 1
        if currentAA != ""
            chunks.push([currentAA, currentMod])
        if currentMod != ""
            mods.push(currentMod)
        if glycan != ""
            glycan = new GlycanComposition(glycan)
        else
            glycan = null
        chunks = (new PeptideSequencePosition(p[0], (m for m in p[1] when m != "")) for p in chunks)
        return [chunks, mods, glycan, nTerm, cTerm]

    constructor: (string) ->
        [@sequence, mods, @glycan, @nTerm, @cTerm] = parser(string)

    format: (colorSource, includeGlycan=true) ->
        positions = []
        if @nTerm == "H"
            nTerm = ""
        else
            nTerm = _formatModification(@nTerm, colorSource).slice(1, -1) + '-'
            positions.push nTerm

        positions = positions.concat (position.format(colorSource) for position in @sequence)

        if @cTerm == "OH"
            cTerm = ""
        else
            cTerm = '-' + (_formatModification(@cTerm, colorSource).slice(1, -1))
            positions.push cTerm
        sequence = positions.join ""
        if includeGlycan
            glycan = @glycan.format(colorSource)
            sequence += ' ' + glycan
        return sequence

