import multiprocessing as mp
from multiprocessing.connection import Connection, wait
from queue import PriorityQueue, Queue
from typing import Callable, List, Optional, Tuple
from enum import Enum


class Command(Enum):
    FINISHED = 1
    READY = 2


class TaskResult:
    def __init__(self, result) -> None:
        self.result = result


class Task:
    def __init__(
        self, func: Callable, args: Tuple, priority: Optional[int] = None
    ) -> None:
        self.func = func
        self.args = args
        self.priority = priority

    def execute(self) -> TaskResult:
        return self.func(*self.args)


class TaskDistributor:
    def __init__(self, num_workers=mp.cpu_count() - 1, use_priority=False) -> None:
        self.task_queue = PriorityQueue() if use_priority else Queue()
        self.num_workers = num_workers
        self.use_priority = use_priority

    def add_task(self, task: Task):
        self.task_queue.put_nowait((task.priority, task))

    def run(self):
        dist_conns: List[Connection] = []
        workers: List[Worker] = []
        for _ in range(self.num_workers):
            dist_conn, work_conn = mp.Pipe()
            dist_conns.append(dist_conn)
            workers.append(Worker(work_conn))

        for worker in workers:
            worker.start()

        available_conns: List[Connection] = []
        busy_conns = list(dist_conns)
        while self.task_queue.qsize() > 0:
            while len(available_conns) > 0 and self.task_queue.qsize() > 0:
                task = self.task_queue.get_nowait()
                conn = available_conns.pop()
                conn.send(task)
                busy_conns.append(conn)

            finished: List = wait(busy_conns)
            for conn in finished:
                result = conn.recv()


class Worker:
    def __init__(self, work_conn: Connection) -> None:
        self.work_conn = work_conn

    def start(self):
        self.process = mp.Process(target=Worker._worker_run, args=(self.work_conn,))

    @staticmethod
    def _worker_run(work_conn: Connection):
        work_conn.send(Command.READY)
        work_conn.poll(timeout=None)
        command = work_conn.recv()
        while command != Command.FINISHED:
            assert isinstance(command, Task)
            result = command.execute()
            work_conn.send(result)
            work_conn.poll(timeout=None)
