"""
Source:

https://lukejharmon.github.io/pcm/chapter10_birthdeath/#section-10.2-the-birth-death-model
"""

from cogent3.core.tree import PhyloNode
from cogent3 import make_tree
import random


def birth_death_tree(
    birth_rate: float, death_rate: float, stopping_time: float, restart_on_fail=False
) -> PhyloNode:
    total_rate = birth_rate + death_rate
    birth_probability = birth_rate / total_rate

    tree = make_tree("(t0:0.0,t1:0.0);")

    tips = list(tree.traverse(self_before=False, self_after=False))

    # print(tips)
    # print([tip.length for tip in tips])

    current_time = 0
    taxa_used = 2
    while current_time < stopping_time:
        waiting_time = random.expovariate(len(tips) * total_rate)
        current_time += waiting_time

        if current_time >= stopping_time:
            adder = stopping_time - (current_time - waiting_time)
            for tip in tips:
                tip.length += adder
            break

        for tip in tips:
            tip.length += waiting_time

        if random.random() < birth_probability:
            # Handle birth event
            birth_tip = tips.pop(random.randrange(0, len(tips)))
            old_tip_name = birth_tip.name
            assert isinstance(old_tip_name, str)
            birth_tip.name = None
            birth_tip.append(make_tree(old_tip_name + ":0.0;"))
            birth_tip.append(make_tree("t" + str(taxa_used) + ":0.0;"))
            taxa_used += 1
            tips.extend(birth_tip.children)
        else:
            # Handle death event
            extinct_tip = tips.pop(random.randrange(0, len(tips)))
            # print("DEATH OF", extinct_tip.name)
            parent: PhyloNode = extinct_tip.parent
            parent.remove_node(extinct_tip)

            assert len(parent.children) == 1
            grandparent: PhyloNode = parent.parent
            if grandparent is None:
                if restart_on_fail:
                    # print("RESTARTING")
                    tree = make_tree("(t0:0.0,t1:0.0);")
                    tips = list(tree.traverse(self_before=False, self_after=False))
                    current_time = 0
                    taxa_used = 2
                    continue
                else:
                    raise RuntimeError("Extinction" + str(birth_probability))
            other_child = parent.children[0]
            other_child.length += parent.length
            if grandparent is None:
                tree = other_child
            else:
                grandparent.remove_node(parent)
                grandparent.append(other_child)

    # print(tree)
    # print([tip.length for tip in tips])

    # print(tree.ascii_art())

    # lengths = []
    # for tip in tree.iter_tips():
    #     lengths.append(
    #         sum(
    #             map(
    #                 lambda x: 0 if x.length is None else x.length,
    #                 tip.ancestors() + [tip],
    #             )
    #         )
    #     )

    # print(lengths)

    return tree


if __name__ == "__main__":
    # random.seed(4)
    tree = birth_death_tree(0.8, 0.2, 1, True)
    print(tree.ascii_art())
