var PeptideSequence, PeptideSequencePosition, _formatModification;

_formatModification = function(modification, colorSource, long) {
  var color, content;
  if (long == null) {
    long = false;
  }
  color = colorSource.get(modification);
  content = escape(long ? modification : modification[0]);
  return "(<span class='modification-chip' style='background-color:" + color + ";padding-left:1px;padding-right:2px;border-radius:2px;' title='" + (escape(modification)) + "' data-modification=" + (escape(modification)) + ">" + content + "</span>)";
};

PeptideSequencePosition = (function() {
  var formatModification;

  formatModification = function(modification, colorSource, long) {
    var color, content;
    if (long == null) {
      long = false;
    }
    color = colorSource.get(modification);
    content = long ? modification : modification[0];
    return "(<span class='modification-chip' style='background-color:" + color + ";padding-left:1px;padding-right:2px;border-radius:2px;' title='" + modification + "' data-modification=" + modification + ">" + content + "</span>)";
  };

  function PeptideSequencePosition(residue, modifications) {
    if (modifications == null) {
      modifications = [];
    }
    this.residue = residue;
    this.modifications = modifications;
  }

  PeptideSequencePosition.prototype.format = function(colorSource) {
    var modifications;
    modifications = this.modifications.map(function(modification) {
      return formatModification(modification, colorSource);
    }).join('');
    return this.residue + modifications;
  };

  return PeptideSequencePosition;

})();

PeptideSequence = (function() {
  var parser;

  parser = function(sequence) {
    var cTerm, chunks, currentAA, currentMod, currentMods, glycan, i, m, mods, nTerm, nextChr, p, parenLevel, state;
    state = "start";
    nTerm = nTerm || "H";
    cTerm = cTerm || "OH";
    mods = [];
    chunks = [];
    glycan = "";
    currentAA = "";
    currentMod = "";
    currentMods = [];
    parenLevel = 0;
    i = 0;
    while (i < sequence.length) {
      nextChr = sequence[i];
      if (nextChr === "(") {
        if (state === "aa") {
          state = "mod";
          parenLevel += 1;
        } else if (state === "start") {
          state = "nTerm";
          parenLevel += 1;
        } else {
          parenLevel += 1;
          if (!((state === "nTerm" || state === "cTerm") && parenLevel === 1)) {
            currentMod += nextChr;
          }
        }
      } else if (nextChr === ")") {
        if (state === "aa") {
          throw new Exception("Invalid Sequence. ) found outside of modification, Position {0}. {1}".format(i, sequence));
        } else {
          parenLevel -= 1;
          if (parenLevel === 0) {
            mods.push(currentMod);
            currentMods.push(currentMod);
            if (state === "mod") {
              state = 'aa';
              if (currentAA === "") {
                chunks.slice(-1)[1] = chunks.slice(-1)[1].concat(currentMods);
              } else {
                chunks.push([currentAA, currentMods]);
              }
            } else if (state === "nTerm") {
              if (sequence[i + 1] !== "-") {
                throw new Exception("Malformed N-terminus for " + sequence);
              }
              nTerm = currentMod;
              state = "aa";
              i += 1;
            } else if (state === "cTerm") {
              cTerm = currentMod;
            }
            currentMods = [];
            currentMod = "";
            currentAA = "";
          } else {
            currentMod += nextChr;
          }
        }
      } else if (nextChr === "|") {
        if (state === "aa") {
          throw new Exception("Invalid Sequence. | found outside of modification");
        } else {
          currentMods.push(currentMod);
          mods.push(currentMod);
          currentMod = "";
        }
      } else if (nextChr === "{") {
        if (state === 'aa' || (state === "cTerm" && parenLevel === 0)) {
          glycan = sequence.slice(i);
          break;
        } else {
          currentMod += nextChr;
        }
      } else if (nextChr === "-") {
        if (state === "aa") {
          state = "cTerm";
          if (currentAA !== "") {
            currentMods.push(currentMod);
            chunks.push([currentAA, currentMods]);
            currentMod = "";
            currentMods = [];
            currentAA = "";
          }
        } else {
          currentMod += nextChr;
        }
      } else if (state === "start") {
        state = "aa";
        currentAA = nextChr;
      } else if (state === "aa") {
        if (currentAA !== "") {
          currentMods.push(currentMod);
          chunks.push([currentAA, currentMods]);
          currentMod = "";
          currentMods = [];
          currentAA = "";
        }
        currentAA = nextChr;
      } else if (state === "nTerm" || state === "mod" || state === "cTerm") {
        currentMod += nextChr;
      } else {
        throw new Exception("Unknown Tokenizer State", currentAA, currentMod, i, nextChr);
      }
      i += 1;
    }
    if (currentAA !== "") {
      chunks.push([currentAA, currentMod]);
    }
    if (currentMod !== "") {
      mods.push(currentMod);
    }
    if (glycan !== "") {
      glycan = new GlycanComposition(glycan);
    } else {
      glycan = null;
    }
    chunks = (function() {
      var j, len, results;
      results = [];
      for (j = 0, len = chunks.length; j < len; j++) {
        p = chunks[j];
        results.push(new PeptideSequencePosition(p[0], (function() {
          var k, len1, ref, results1;
          ref = p[1];
          results1 = [];
          for (k = 0, len1 = ref.length; k < len1; k++) {
            m = ref[k];
            if (m !== "") {
              results1.push(m);
            }
          }
          return results1;
        })()));
      }
      return results;
    })();
    return [chunks, mods, glycan, nTerm, cTerm];
  };

  function PeptideSequence(string) {
    var mods, ref;
    ref = parser(string), this.sequence = ref[0], mods = ref[1], this.glycan = ref[2], this.nTerm = ref[3], this.cTerm = ref[4];
  }

  PeptideSequence.prototype.format = function(colorSource, includeGlycan) {
    var cTerm, glycan, nTerm, position, positions, sequence;
    if (includeGlycan == null) {
      includeGlycan = true;
    }
    positions = [];
    if (this.nTerm === "H") {
      nTerm = "";
    } else {
      nTerm = _formatModification(this.nTerm, colorSource).slice(1, -1) + '-';
      positions.push(nTerm);
    }
    positions = positions.concat((function() {
      var j, len, ref, results;
      ref = this.sequence;
      results = [];
      for (j = 0, len = ref.length; j < len; j++) {
        position = ref[j];
        results.push(position.format(colorSource));
      }
      return results;
    }).call(this));
    if (this.cTerm === "OH") {
      cTerm = "";
    } else {
      cTerm = '-' + (_formatModification(this.cTerm, colorSource).slice(1, -1));
      positions.push(cTerm);
    }
    sequence = positions.join("");
    if (includeGlycan) {
      glycan = this.glycan.format(colorSource);
      sequence += ' ' + glycan;
    }
    return sequence;
  };

  return PeptideSequence;

})();

//# sourceMappingURL=peptide-sequence.js.map
