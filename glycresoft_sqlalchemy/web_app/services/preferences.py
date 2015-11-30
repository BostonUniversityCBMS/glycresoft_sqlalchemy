from flask import Blueprint, g, jsonify, render_template, request


app_config = Blueprint("preferences", __name__)


default_preferences = {
    "minimum_ms2_score": 0.4
}


@app_config.route("/preferences")
def show_preferences():
    preferences = g.manager.app_data.get("preferences", default_preferences)
    for k, v in default_preferences.items():
        if k not in preferences:
            preferences[k] = v
    return render_template("components/preferences.templ", **preferences)


@app_config.route("/preferences", methods=["POST"])
def update_preferences():
    new_preferences = request.values.to_dict()
    print "new_preferences"
    try:
        preferences = g.manager["preferences"]
    except KeyError:
        preferences = {}
    for k, v in default_preferences.items():
        if k not in preferences:
            preferences[k] = v
    preferences.update(new_preferences)
    g.manager["preferences"] = preferences
    print("\n%r\n" % g.manager["preferences"])
    return jsonify(**dict(preferences.items()))
