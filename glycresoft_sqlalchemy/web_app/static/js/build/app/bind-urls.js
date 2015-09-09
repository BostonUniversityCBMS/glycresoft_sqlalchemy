var ActionBook, DataSource, PartialSource, makeAPIGet, makePartialGet;

ActionBook = {
  home: {
    container: '#home-layer',
    name: 'home-layer'
  },
  addSample: {
    contentURL: '/add_sample',
    name: 'add-sample'
  },
  peakGroupingMatchSamples: {
    contentURL: '/peak_grouping_match_samples',
    name: "peak-grouping-match-samples"
  },
  tandemMatchSamples: {
    contentURL: '/tandem_match_samples',
    name: 'tandem-match-samples'
  },
  naiveGlycopeptideSearchSpace: {
    contentURL: "/glycopeptide_search_space",
    name: "glycopeptide-search-space"
  },
  naiveGlycanSearchSpace: {
    contentURL: "/glycan_search_space",
    name: "glycan-search-space"
  },
  viewDatabaseSearchResults: {
    contentURLTemplate: "/view_database_search_results/{hypothesis_sample_match_id}",
    name: "view-database-search-results",
    method: "post"
  }
};

makeAPIGet = function(url) {
  return function(callback) {
    return $.get(url).success(callback);
  };
};

DataSource = {
  hypotheses: makeAPIGet("/api/hypotheses"),
  samples: makeAPIGet("/api/samples"),
  hypothesisSampleMatches: makeAPIGet("/api/hypothesis_sample_matches"),
  tasks: makeAPIGet("/api/tasks")
};

makePartialGet = function(url, method) {
  return function(parameters, callback) {
    return $[method](url.format(parameters)).success(callback);
  };
};

PartialSource = {
  glycopeptideCompositionDetailsModal: makePartialGet('/view_database_search_results/view_glycopeptide_composition_details/{id}', "get")
};

//# sourceMappingURL=bind-urls.js.map
