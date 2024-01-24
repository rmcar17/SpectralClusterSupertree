# print("beginning import")
import multiprocessing as mp
from multiprocessing.connection import Connection, wait
from queue import PriorityQueue, Queue
from cogent3 import TreeNode
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple, Union
from enum import Enum
import numpy as np
from spectral_cluster_supertree.scs.scs import (
    _component_to_names_set,
    _connect_trees,
    _contract_proper_cluster_graph,
    _denamify,
    _generate_induced_trees_with_weights,
    _get_all_tip_names,
    _get_graph_components,
    _proper_cluster_graph_edges,
    _tip_names_to_tree,
    spectral_cluster_graph,
    spectral_cluster_supertree,
)
from cogent3.core.tree import TreeBuilder

# print("end import")

taskID = int


def parallel_scs_split(
    trees: Sequence[TreeNode],
    pcg_weighting: str = "one",
    normalise_pcg_weights: bool = False,
    depth_normalisation: bool = False,
    contract_edges: bool = True,
    weights: Optional[Sequence[float]] = None,
) -> Union[TreeNode, "MergeResultOfTasks"]:
    """
    Spectral Cluster Supertree (SCS).

    Constructs a supertree from a collection of input trees. The supertree
    method is inspired by Min-Cut Supertree (Semple & Steel, 2000), using
    spectral clustering instead of min-cut to improve efficiency.

    The set of input trees must overlap, the optional weights parameter
    allows the biasing of some trees over others.

    Args:
        trees (Sequence[TreeNode]): Overlapping subtrees.
        weights (Optional[Sequence[float]]): Optional weights for the trees.

    Returns:
        TreeNode: The supertree containing all taxa in the input trees.
    """

    assert len(trees) >= 1, "there must be at least one tree"

    assert pcg_weighting in ["one", "branch", "depth"]

    # Input trees are of equal weight if none is specified
    if weights is None:
        weights = [1.0 for _ in range(len(trees))]

    assert len(trees) == len(weights), "trees and weights must be of same length"

    if len(trees) == 1:  # If there is only one tree left, we can simply graft it on
        _denamify(trees[0])
        return trees[0]

    # The vertices of the proper cluster graph
    # are the names of the tips of all trees
    all_names = _get_all_tip_names(trees)

    if len(all_names) < 100:
        return spectral_cluster_supertree(
            trees,
            pcg_weighting=pcg_weighting,
            normalise_pcg_weights=normalise_pcg_weights,
            depth_normalisation=depth_normalisation,
            contract_edges=contract_edges,
            weights=weights,
        )

    # If there are less than or only two names, can instantly return a tree
    if len(all_names) <= 2:
        tree = _tip_names_to_tree(all_names)
        return tree

    pcg_vertices = set((name,) for name in all_names)

    (
        pcg_edges,
        pcg_weights,
        taxa_ocurrences,
        taxa_co_occurrences,
    ) = _proper_cluster_graph_edges(
        pcg_vertices,
        trees,
        weights,
        pcg_weighting,
        normalise_pcg_weights,
        depth_normalisation,
    )

    components = _get_graph_components(pcg_vertices, pcg_edges)

    if len(components) == 1:
        if contract_edges:
            # Modifies the proper cluster graph inplace
            _contract_proper_cluster_graph(
                pcg_vertices,
                pcg_edges,
                pcg_weights,
                taxa_ocurrences,
                taxa_co_occurrences,
            )
        components = spectral_cluster_graph(
            pcg_vertices, pcg_weights, np.random.RandomState()
        )

    return MergeResultOfTasks(
        [
            Task(
                "handle_component",
                (
                    component,
                    trees,
                    weights,
                    pcg_weighting,
                    normalise_pcg_weights,
                    depth_normalisation,
                    contract_edges,
                ),
            )
            for component in components
        ],
        "_connect_trees_fix_tips",
        kwargs={"all_tip_names": all_names},
    )


def handle_component(
    component,
    trees,
    weights,
    pcg_weighting,
    normalise_pcg_weights,
    depth_normalisation,
    contract_edges,
) -> Union[TreeNode, "MergeResultOfTasks"]:
    component = _component_to_names_set(component)
    # Trivial case for if the size of the component is <=2
    # Simply add a tree expressing that
    if len(component) <= 2:
        return _tip_names_to_tree(component)

    # Otherwise, need to induce the trees on each compoment
    # and recursively call SCS

    # Note, inducing could possible remove trees.
    new_induced_trees, new_weights = _generate_induced_trees_with_weights(
        component, trees, weights
    )

    # Find the supertree for the induced trees
    child_tree = parallel_scs_split(
        new_induced_trees,
        pcg_weighting,
        normalise_pcg_weights,
        depth_normalisation,
        contract_edges,
        new_weights,
    )

    return child_tree


def _connect_trees_fix_tips(*args, **kwargs) -> TreeNode:
    """
    Connects the input trees by making them adjacent to a new root.

    Args:
        trees (Iterable[TreeNode]): The input trees to connect

    Returns:
        TreeNode: A tree connecting all the input trees
    """
    trees = list(args)

    all_tip_names = kwargs["all_tip_names"]

    for tree in trees:
        all_tip_names.difference_update(tree.get_tip_names())

    trees.extend(map(lambda x: _tip_names_to_tree((x,)), all_tip_names))

    if len(trees) == 1:
        (one,) = trees  # Unpack only tree
        return one
    tree_builder = TreeBuilder(constructor=TreeNode).edge_from_edge  # type: ignore
    return tree_builder(None, trees)


class Command(Enum):
    FINISHED = 1
    READY = 2


class TaskResult:
    def __init__(self, task_id: taskID, result) -> None:
        self.task_id = task_id
        self.result = result


class Task:
    number_of_tasks: taskID = 0

    funcs = {
        "parallel_scs_split": parallel_scs_split,
        "handle_component": handle_component,
        "_connect_trees_fix_tips": _connect_trees_fix_tips,
    }

    def __init__(
        self,
        func: str,
        args: Optional[Sequence],
        priority: Optional[int] = None,
        kwargs=None,
    ) -> None:
        self.func = Task.funcs[func]
        self.args = args
        self.priority = priority
        self.dependent_tasks = set()
        self.task_id = None
        self.kwargs = kwargs or {}

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
        return TaskResult(self.task_id, self.func(*self.args, **self.kwargs))

    def __str__(self) -> str:
        return (
            "Call "
            + str(self.func)
            + " on "
            + str(self.args)
            + f" (tid: {self.task_id})"
        )


class MergeResultOfTasks:
    def __init__(self, subtasks: List[Task], merger: str, kwargs=None) -> None:
        self.subtasks = subtasks
        self.merger = merger
        self.subtask_ids: Optional[List[taskID]] = None
        self.kwargs = kwargs or {}

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
        task = Task(self.merger, None, kwargs=self.kwargs)
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
        self.process = mp.Process(
            target=Worker._worker_run, args=(self.work_conn,), daemon=False
        )

    def start(self):
        # print("STARTING PROCESS")
        self.process.start()
        # print("STARTED")

    @staticmethod
    def _worker_run(work_conn: Connection):
        # print("WORKER SENDING READY", __name__)
        from cogent3 import make_tree

        work_conn.send(Command.READY)
        # print("WAITING FOR COMMAND")
        work_conn.poll(timeout=None)
        command = work_conn.recv()
        while command != Command.FINISHED:
            # print("COMMAND", command)
            assert isinstance(command, Task)
            result = command.execute()

            work_conn.send(result)
            work_conn.poll(timeout=None)
            command = work_conn.recv()


# print("loaded")
