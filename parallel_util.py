import multiprocessing as mp
from multiprocessing.connection import Connection, wait
from queue import PriorityQueue, Queue
import time
from typing import Callable, Dict, List, Optional, Set, Tuple
from enum import Enum

taskID = int


class Command(Enum):
    FINISHED = 1
    READY = 2


class TaskResult:
    def __init__(self, task_id: taskID, result) -> None:
        self.task_id = task_id
        self.result = result


class Task:
    number_of_tasks: taskID = 0

    def __init__(
        self, func: Callable, args: Optional[Tuple], priority: Optional[int] = None
    ) -> None:
        self.func = func
        self.args = args
        self.priority = priority
        self.dependent_tasks = set()
        self.task_id = None

    def distributor_add_dependent_task(self, task: taskID):
        self.dependent_tasks.add(task)

    def distributor_set_task_id(self):
        self.task_id = Task.number_of_tasks
        Task.number_of_tasks += 1

    def execute(self) -> TaskResult:
        assert self.task_id is not None
        assert self.args is not None
        return TaskResult(self.task_id, self.func(*self.args))

    def __str__(self) -> str:
        return "Call " + str(self.func) + " on " + str(self.args)


class MergeResultOfTasks:
    def __init__(self, subtasks: List[Task], merger: Callable) -> None:
        self.subtasks = subtasks
        self.merger = merger
        self.subtask_ids: Optional[List[taskID]] = None

    def get_subtask_ids(self) -> List[taskID]:
        assert self.subtask_ids is not None
        return self.subtask_ids

    def distributor_set_task_ids(self):
        self.subtask_ids = []
        for subtask in self.subtasks:
            subtask.distributor_set_task_id()
            assert subtask.task_id is not None
            self.subtask_ids.append(subtask.task_id)

    def distributor_generate_dependent_task(self) -> Task:
        task = Task(self.merger, None)
        task.distributor_set_task_id()
        assert task.task_id is not None
        for subtask in self.subtasks:
            subtask.distributor_add_dependent_task(task.task_id)
        return task


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

        saved_results: Dict[taskID, TaskResult] = {}
        dependencies: Dict[taskID, Set[taskID]] = {}
        dependent_tasks: Dict[taskID, Task] = {}

        available_conns: List[Connection] = []
        busy_conns = list(dist_conns)
        result = None
        while self.task_queue.qsize() > 0 or len(busy_conns) > 0:
            print(len(available_conns), self.task_queue.qsize())
            while len(available_conns) > 0 and self.task_queue.qsize() > 0:
                task: Task
                priority, task = self.task_queue.get_nowait()
                if task.task_id is None:
                    task.distributor_set_task_id()
                conn = available_conns.pop()
                conn.send(task)
                busy_conns.append(conn)

            finished = wait(busy_conns)
            for conn in finished:
                assert isinstance(conn, Connection)
                result = conn.recv()

                if isinstance(result, MergeResultOfTasks):
                    result.distributor_set_task_ids()
                    dependent_task = result.distributor_generate_dependent_task()
                    dependent_task_id = dependent_task.task_id
                    assert dependent_task_id is not None

                    assert dependent_task_id not in dependent_tasks
                    dependent_tasks[dependent_task_id] = dependent_task

                    assert dependent_task_id not in dependencies
                    dependencies[dependent_task_id] = set(result.get_subtask_ids())

                    for task in result.subtasks:
                        assert (
                            task.task_id not in dependent_tasks
                            and task.task_id is not None
                        )
                        dependent_tasks[task.task_id] = task
                        self.add_task(task)
                elif isinstance(result, TaskResult):
                    if result.task_id in dependent_tasks:
                        relevant_task = dependent_tasks.pop(result.task_id)
                        for dependent_id in relevant_task.dependent_tasks:
                            dependencies
                            # TODO FROM HERE
                busy_conns.remove(conn)
                available_conns.append(conn)

        for conn in dist_conns:
            conn.send(Command.FINISHED)

        return result


class Worker:
    def __init__(self, work_conn: Connection) -> None:
        self.work_conn = work_conn
        self.process = mp.Process(target=Worker._worker_run, args=(self.work_conn,))

    def start(self):
        self.process.start()

    @staticmethod
    def _worker_run(work_conn: Connection):
        work_conn.send(Command.READY)
        work_conn.poll(timeout=None)
        command = work_conn.recv()
        while command != Command.FINISHED:
            print(command)
            assert isinstance(command, Task)
            result = command.execute()

            work_conn.send(result)
            start = time.time()
            work_conn.poll(timeout=None)
            command = work_conn.recv()
