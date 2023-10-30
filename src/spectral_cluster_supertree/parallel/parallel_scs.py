from multiprocessing.connection import Connection, wait
from typing import List, Optional, Sequence, Union
from cogent3.core.tree import TreeNode
import multiprocessing as mp
from parallel_util import MergeResultOfTasks, Task, TaskDistributor, TaskResult
from spectral_cluster_supertree import scs
from cogent3 import make_tree


def func(x, y):
    return x + y + 1


def sum_sequence(sequence: Sequence[int]) -> int:
    if len(sequence) == 1:
        return sequence[0]
    midpoint = len(sequence) // 2
    return sum_sequence(sequence[:midpoint]) + sum_sequence(sequence[midpoint:])


def parallel_sum_sequence_split(
    sequence: Sequence[int],
) -> Union[int, MergeResultOfTasks]:
    if len(sequence) == 1:
        return sequence[0]
    midpoint = len(sequence) // 2
    return MergeResultOfTasks(
        [
            Task(parallel_sum_sequence_split, (sequence[:midpoint],)),
            Task(parallel_sum_sequence_split, (sequence[midpoint:],)),
        ],
        parallel_sum_sequence_add,
    )


def parallel_sum_sequence_add(left: int, right: int) -> int:
    return left + right


def parallel_spectral_cluster_supertree(
    trees: Sequence[TreeNode],
    weights: Optional[Sequence[float]] = None,
    num_workers: int = mp.cpu_count() - 1,
) -> TreeNode:
    if num_workers == 0:
        return scs(trees, weights)

    print("INIT DIST")
    distributor = TaskDistributor(num_workers)
    # distributor.add_task(Task(func, (3, 4)))
    distributor.initialise_workers()
    n = 10000
    seq = list(range(n + 1))
    distributor.add_task(Task(parallel_sum_sequence_split, (seq,)))
    import time

    start = time.time()
    print(distributor.run(), time.time() - start)
    start = time.time()
    print(my_sum(seq), time.time() - start)
    print((n * (n + 1)) // 2)


def my_sum(seq):
    total = 0
    for s in seq:
        total += s
    return total


if __name__ == "__main__":
    parallel_spectral_cluster_supertree(
        [make_tree("(a,(b,c))"), make_tree("(b,(a,d))")], num_workers=11
    )
