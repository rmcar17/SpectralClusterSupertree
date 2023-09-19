import multiprocessing as mp
from multiprocessing.connection import Connection, wait
from queue import PriorityQueue, Queue
import time
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple
from enum import Enum

taskID = int
# print("STARTING")


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
        self, func: Callable, args: Optional[Sequence], priority: Optional[int] = None
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

    def add_argument(self, arg):
        # Assumes merge function is unordered in parameters. May want to make more general
        if self.args is None:
            self.args = list()
        if not isinstance(self.args, list):
            self.args = list(self.args)
        self.args.append(arg)

    def execute(self) -> TaskResult:
        assert self.task_id is not None
        assert self.args is not None
        # print(self.func, self.args)
        return TaskResult(self.task_id, self.func(*self.args))

    def __str__(self) -> str:
        return (
            "Call "
            + str(self.func)
            + " on "
            + str(self.args)
            + f" (tid: {self.task_id})"
        )


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
        self.workers_initialised = False

    def add_task(self, task: Task):
        self.task_queue.put_nowait((task.priority, task))

    def initialise_workers(self):
        self.dist_conns: List[Connection] = []
        self.workers: List[Worker] = []
        for _ in range(self.num_workers):
            dist_conn, work_conn = mp.Pipe()
            self.dist_conns.append(dist_conn)
            self.workers.append(Worker(work_conn))

        for worker in self.workers:
            worker.start()
        self.workers_initialised = True

    def run(self):
        assert self.workers_initialised

        dependencies: Dict[taskID, Set[taskID]] = {}
        dependent_tasks: Dict[taskID, Task] = {}

        available_conns: List[Connection] = []
        busy_conns = list(self.dist_conns)
        result = None
        last_task_result = None
        while self.task_queue.qsize() > 0 or len(busy_conns) > 0:
            # print(len(available_conns), self.task_queue.qsize())
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

                busy_conns.remove(conn)
                available_conns.append(conn)

                if result == Command.READY:
                    continue

                assert isinstance(result, TaskResult)

                merge_result = result.result

                if isinstance(merge_result, MergeResultOfTasks):
                    # print(
                    #     "GOT MERGE RESULT",
                    #     result.task_id,
                    #     None
                    #     if result.task_id not in dependent_tasks
                    #     else dependent_tasks[result.task_id].dependent_tasks,
                    # )
                    merge_result.distributor_set_task_ids()
                    dependent_task = merge_result.distributor_generate_dependent_task()
                    dependent_task_id = dependent_task.task_id
                    assert dependent_task_id is not None

                    assert dependent_task_id not in dependent_tasks
                    dependent_tasks[dependent_task_id] = dependent_task

                    assert dependent_task_id not in dependencies
                    dependencies[dependent_task_id] = set(
                        merge_result.get_subtask_ids()
                    )

                    assert result.task_id is not None
                    # print("HANDLING dependencies")
                    if result.task_id in dependent_tasks:
                        finished_task = dependent_tasks.pop(
                            result.task_id
                        )  # The task the task result is from
                        # print(finished_task)
                        # print(finished_task.dependent_tasks)
                        # print(merge_result.subtasks)

                        # For each task which has a dependency on the finished task
                        for dependent_id in finished_task.dependent_tasks:
                            dependencies[dependent_id].remove(result.task_id)
                            # o_dependent_task = dependent_tasks[dependent_id]
                            # o_dependent_task.add_argument(result.result)
                            assert dependent_task.task_id is not None
                            dependencies[dependent_id].add(dependent_task.task_id)
                            dependent_tasks[dependent_task.task_id].dependent_tasks.add(
                                dependent_id
                            )

                    for task in merge_result.subtasks:
                        assert (
                            task.task_id not in dependent_tasks
                            and task.task_id is not None
                        )
                        dependent_tasks[task.task_id] = task
                        self.add_task(task)
                else:
                    last_task_result = result
                    # print(dependencies)
                    if result.task_id in dependent_tasks:
                        finished_task = dependent_tasks.pop(
                            result.task_id
                        )  # The task the task result is from

                        # For each task which has a dependency on the finished task
                        for dependent_id in finished_task.dependent_tasks:
                            dependencies[dependent_id].remove(result.task_id)
                            dependent_task = dependent_tasks[dependent_id]
                            dependent_task.add_argument(result.result)

                            # The task is ready to be dispatched
                            if len(dependencies[dependent_id]) == 0:
                                del dependencies[dependent_id]
                                self.add_task(dependent_task)

        for conn in self.dist_conns:
            conn.send(Command.FINISHED)

        # assert isinstance(result, TaskResult)
        # print(result, last_task_result)
        if last_task_result is None:
            return None
        # print(last_task_result, last_task_result.result)
        return last_task_result.result


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
            # print(command)
            assert isinstance(command, Task)
            result = command.execute()

            work_conn.send(result)
            work_conn.poll(timeout=None)
            command = work_conn.recv()
