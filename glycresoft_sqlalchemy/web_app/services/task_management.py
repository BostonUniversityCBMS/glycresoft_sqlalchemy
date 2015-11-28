from flask import Response, g, Blueprint

task_actions = Blueprint("task_management", __name__)


@task_actions.route("/internal/log/<task_id>")
def send_log(task_id):
    return Response("<pre>%s</pre>" % open(
        g.manager.get_task_path(task_id + '.log'), 'r').read(), mimetype='application/text')
