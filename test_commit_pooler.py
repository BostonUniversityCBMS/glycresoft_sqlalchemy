from glycresoft_sqlalchemy.data_model import commit_pool as commit_pooler
from multiprocessing import Process
import time

if __name__ == '__main__':
    commit_pooler.DatabaseManager("./datafiles/test.db", clear=True).initialize()
    pooler = commit_pooler.CommitPooler("./datafiles/test.db")
    task_queue = pooler.start()
    print "Pooler started"
    commit_pooler.tester(task_queue)
    print "Sent Data"
    for i in range(6):
        proc = Process(target=commit_pooler.tester, args=(task_queue,))
        proc.start()
    time.sleep(5)
    pooler.terminate()
    print "Checking Database"
    s = commit_pooler.DatabaseManager(pooler.database_path).session()
    print(list(s.query(commit_pooler.Experiment)))
