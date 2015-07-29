ActionBook =
    home:
        container: '#home-layer'
        name: 'home-layer'
    addSample:
        contentURL: '/add_sample'
        name: 'add-sample'
    tandemMatchSamples:
        contentURL: '/tandem_match_samples'
        name: 'tandem-match-samples'
    naiveGlycopeptideSearchSpace:
        contentURL: "/glycopeptide_search_space"
        name: "glycopeptide-search-space"
    naiveGlycanSearchSpace:
        contentURL: "/glycan_search_space"
        name: "glycan-search-space"
    viewDatabaseSearchResults:
        contentURLTemplate: "/view_database_search_results/{}"
        name: "view-database-search-results"

makeAPIGet = (url) -> (callback) -> $.get(url).success(callback)

DataSource =
    hypotheses: makeAPIGet "/api/hypotheses"
    samples: makeAPIGet "/api/samples"
    hypothesisSampleMatches: makeAPIGet "/api/hypothesis_sample_matches"
    tasks: makeAPIGet "/api/tasks"
