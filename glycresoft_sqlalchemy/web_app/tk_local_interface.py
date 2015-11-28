import os
import Tkinter as tk
import tkFileDialog
import tkMessageBox

from glycresoft_sqlalchemy.utils import sqlitedict, appdir
from glycresoft_sqlalchemy.web_app.serve import run

print "Imported"
root = tk.Tk()
root.withdraw()


def get_directory(**options):
    return tkFileDialog.askdirectory(**options)


def get_open_file(**options):
    return tkFileDialog.askopenfilename(**options)


def get_save_file(**options):
    return tkFileDialog.asksaveasfilename(**options)


ask_yesnocancel = tkMessageBox.askyesnocancel
show_error = tkMessageBox.showerror
show_info = tkMessageBox.showinfo
show_warning = tkMessageBox.showwarning


DEFAULT_PORT = 5953
DEFAULT_EXTERNAL = False
DEFAULT_NO_EXECUTE_TASKS = False


class ApplicationStartup(object):

    def __init__(self, history_file=None):
        self.history_file = (history_file)
        self.history = sqlitedict.open(history_file)
        self.project_paths = []
        self.last_project = None
        self.project = None
        self.port = DEFAULT_PORT
        self.no_execute_tasks = DEFAULT_NO_EXECUTE_TASKS
        self.external = DEFAULT_EXTERNAL

    def save(self):
        self.history["project_paths"] = self.project_paths
        self.history["last_project"] = self.last_project

        self.history["port"] = self.port
        self.history["no_execute_tasks"] = self.no_execute_tasks
        self.history["external"] = self.external
        self.history.commit()

    def get_project_path(self):
        pass

    def get_new_project_path(self):
        new_project_path = get_directory(title="Create a new project", initialdir=appdir.user_data_dir())
        if new_project_path == "":
            new_project_path = None
        return new_project_path

    def load_history(self):
        self.project_paths = self.history.get("project_paths", [])
        self.last_project = self.history.get("last_project")
        self.port = self.history.get("port", DEFAULT_PORT)
        self.no_execute_tasks = self.history.get("no_execute_tasks", DEFAULT_NO_EXECUTE_TASKS)
        self.external = self.history.get("external", DEFAULT_EXTERNAL)

        self.last_project = new_project_path = self.get_new_project_path()
        self.project_paths.append(new_project_path)
        if new_project_path is None:
            raise ValueError("New Project Path Is None")
        self.save()
        self.project = Project(self.last_project)


def startup():
    app = ApplicationStartup()
    app.load_history()
    run(str(app.project.store_file), app.external, app.no_execute_tasks, app.port)


class Project(object):
    main_project_file_name = "store.glycresoftdb"

    def __init__(self, directory):
        self.directory = directory
        self.store_file = os.path.join(self.directory, self.main_project_file_name)

    def __repr__(self):
        return "Project(directory=%s)" % self.directory

    def __str__(self):
        return self.store_file


if __name__ == '__main__':
    startup()
