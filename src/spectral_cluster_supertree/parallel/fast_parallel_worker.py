# import multiprocessing as mp
# from multiprocessing.connection import Connection, wait
# from enum import Enum
# from typing import Optional, Sequence
# from cogent3 import TreeNode


# class Command(Enum):
#     FINISHED = 1
#     READY = 2


# class Job:
#     pass


# class SplitJob(Job):
#     def __init__(
#         self,
#         trees: Sequence[TreeNode],
#         pcg_weighting: str,
#         normalise_pcg_weights: bool,
#         depth_normalisation: bool,
#         contract_edges: bool,
#         weights: Optional[Sequence[float]],
#     ) -> None:
#         super().__init__()
#         self.trees = trees
#         self.pcg_weighting = pcg_weighting
#         self.normalise_pcg_weights = normalise_pcg_weights
#         self.depth_normalisation = depth_normalisation
#         self.contract_edges = contract_edges
#         self.weights = weights


# class MergeJob(Job):
#     def __init__(self) -> None:
#         super().__init__()


# def worker(conn: Connection):
#     conn.send(Command.READY)
#     conn.poll(timeout=None)

#     command = conn.recv()
#     while command != Command.FINISHED:
#         assert isinstance(command, Job)

#         if isinstance(command, SplitJob):
#             result = scs_split(command)
#         else:
#             isinstance(
#                 command,
#             )
