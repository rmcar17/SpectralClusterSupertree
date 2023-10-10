import os
import subprocess
from typing import List, Optional, Tuple
from cogent3 import make_tree
from cogent3.core.tree import TreeNode
import time
from distance import (
    cluster_matching_distance,
    grf_distance,
    grf_distance_slow,
    rf_distance,
)

data_path = "data/superfine/"
model_trees_end = ".model_tree"
source_trees_end = ".source_trees"


# def simulated_experiment(taxa, density):
#     path = data_path + f"{taxa}-taxa/{density}/"
#     data_file_starts = []
#     for file in os.listdir(path):
#         if file.endswith(source_trees_end):
#             data_file_starts.append(".".join(file.split(".")[:-1]))

#     for data_file_start in data_file_starts:
#         source_tree_file = path + data_file_start + source_trees_end
#         model_tree_file = path + data_file_start + model_trees_end
#         result_sup = subprocess.run(
#             ["./run_sup.sh", source_tree_file],
#             stdout=subprocess.PIPE,
#             stderr=subprocess.PIPE,
#         )

#         result_scs = subprocess.run(
#             ["./run_scs.sh", source_tree_file],
#             stdout=subprocess.PIPE,
#             stderr=subprocess.PIPE,
#         )
#         # result = subprocess.run(["ls", "-l"], stdout=subprocess.PIPE)
#         # print(result.stdout)
#         # print(result.stderr)
#         sup_tree = (
#             make_tree(result_sup.stdout.decode("utf-8").strip())
#             .bifurcating()
#             .unrooted()
#         )
#         scs_tree = (
#             make_tree(result_scs.stdout.decode("utf-8").strip())
#             .bifurcating()
#             .unrooted()
#         )
#         print("SUPERFINE:", sup_tree)
#         print("SCS:", scs_tree)

#         with open(model_tree_file, "r") as f:
#             model = make_tree(f.read().strip()).bifurcating().unrooted()

#         print(sup_tree.lin_rajan_moret(scs_tree))
#         print(model.lin_rajan_moret(sup_tree))
#         print(model.lin_rajan_moret(scs_tree))

#         break

SUP = "SUP"
SCS = "SCS"
MCS = "MCS"
BCD = "BCD"

SCRIPTS = {
    SUP: "./run_sup.sh",
    SCS: "./run_scs.sh",
    MCS: "./run_mcs.sh",
    BCD: "./run_bcd.sh",
}


def run_methods(
    source_tree_file: str, model_tree_file: str, methods: List[str], verbosity=1
):
    results = {}
    for method in methods:
        if verbosity >= 2:
            print("Running Method", method)
        start_time = time.time()
        result = subprocess.run(
            [SCRIPTS[method], source_tree_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        end_time = time.time()
        try:
            tree = make_tree(result.stdout.decode("utf-8").strip()).bifurcating()
        except:
            tree = None
        results[method] = (tree, end_time - start_time)

    with open(model_tree_file, "r") as f:
        model = make_tree(f.read().strip()).bifurcating()

    for method in methods:
        tree, time_result = results[method]
        if verbosity >= 1:
            if verbosity >= 2:
                print("Computing distances for", method)
            if tree is None:
                print(f"{method}: time={time_result:.2f}s failed: {source_tree_file}) ")
            else:
                rf = rf_distance(model, tree)
                grf = grf_distance(model, tree)
                mat = cluster_matching_distance(model, tree)
                print(f"{method}: time={time_result:.2f}s RF={rf} MAT={mat} GRF={grf}")

    return results


def report(source_tree_file, model_tree_file, min_cut=False, verbose=False):
    # sup_start = time.time()
    # result_sup = subprocess.run(
    #     ["./run_sup.sh", source_tree_file],
    #     stdout=subprocess.PIPE,
    #     stderr=subprocess.PIPE,
    # )
    # sup_time = time.time() - sup_start

    scs_start = time.time()
    result_scs = subprocess.run(
        ["./run_scs.sh", source_tree_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    scs_time = time.time() - scs_start

    bcd_start = time.time()
    result_bcd = subprocess.run(
        ["./run_bcd.sh", source_tree_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    # print(result_bcd.stderr)
    # print(result_bcd.stdout)
    bcd_time = time.time() - bcd_start

    mcs_tree = None
    mcs_time = None
    if min_cut:
        mcs_start = time.time()
        result_mcs = subprocess.run(
            ["./run_mcs.sh", source_tree_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        mcs_time = time.time() - mcs_start
        mcs_tree = (
            make_tree(result_mcs.stdout.decode("utf-8").strip()).bifurcating()
            # .unrooted()
        )

    # try:
    #     sup_tree = (
    #         make_tree(result_sup.stdout.decode("utf-8").strip()).bifurcating()
    #         # .unrooted()
    #     )
    # except:
    #     sup_tree = None
    #     print(result_sup.stdout.decode("utf-8").strip())
    #     print(result_sup.stderr.decode("utf-8").strip())

    scs_tree = make_tree(
        result_scs.stdout.decode("utf-8").strip()
    ).bifurcating()  # .unrooted()

    bcd_tree = make_tree(result_bcd.stdout.decode("utf-8").strip()).bifurcating()

    if verbose:
        # print("SUPERFINE:", sup_tree)
        print("SCS:", scs_tree)
        print("BCD", bcd_tree)
        if min_cut:
            print("MCS:", mcs_tree)

    with open(model_tree_file, "r") as f:
        model = make_tree(f.read().strip()).bifurcating()  # .unrooted()

    # if sup_tree is not None:
    #     # print(sup_tree.lin_rajan_moret(scs_tree))
    #     print(
    #         f"SUP: time={sup_time:.2f} RF={rf_distance(model, sup_tree)} GRF={grf_distance(model, sup_tree)}, {grf_distance_slow(model, sup_tree)}, MAT={cluster_matching_distance(model, sup_tree)}"
    #     )
    print(
        f"SCS: time={scs_time:.2f} RF={rf_distance(model, scs_tree)} GRF={grf_distance(model, scs_tree)}, {grf_distance_slow(model, scs_tree)}, MAT={cluster_matching_distance(model, scs_tree)}"
    )
    print(
        f"BCD: time={bcd_time:.2f} RF={rf_distance(model, bcd_tree)} GRF={grf_distance(model, bcd_tree)}, {grf_distance_slow(model, bcd_tree)}, MAT={cluster_matching_distance(model, bcd_tree)}"
    )
    # try:
    # except AssertionError as e:
    #     print(
    #         f"BCD generated tree with {len(set(model.get_tip_names()).difference(bcd_tree.get_tip_names()))} missing tips"
    #     )
    #     print(e)
    if mcs_tree is not None and mcs_time is not None:
        print(
            f"MCS: time={mcs_time:.2f} RF={rf_distance(model, mcs_tree)} GRF={grf_distance(model, mcs_tree)}, {grf_distance_slow(model, mcs_tree)}, MAT={cluster_matching_distance(model, mcs_tree)}"
        )


def run_experiment_super_triplets(d: int, k: int, methods: List):
    assert d in [25, 50, 75]
    assert k in [10, 20, 30, 40, 50]

    number_of_experiments = 100

    for i in range(number_of_experiments):
        tree_number = i + 1
        source_file = f"bcd/data/SuperTripletsBenchmark/source-trees/d{d}/k{k}/data-d{d}-k{k}-{tree_number}_phybp-s0r.nwk.source_trees"
        model_file = f"bcd/data/SuperTripletsBenchmark/model-trees/model-{tree_number}.nwk.model_tree"
        print(f"Results for {i} ({source_file}):")
        run_methods(source_file, model_file, methods)


def run_experiment_smidgen(taxa: int, density: int, methods: List):
    assert taxa in [100, 500, 1000]
    assert density in [20, 50, 75, 100]

    number_of_experiments = 10 if taxa == 1000 else 30

    for i in range(number_of_experiments):
        file = f"data/superfine/{taxa}-taxa/{density}/sm_data.{i}"
        print(f"Results for {i} ({file}):")
        run_methods(file + ".source_trees", file + ".model_tree", methods)


def run_experiment_smidgen_og(
    taxa: int, density: int, methods: List, verbosity: int = 1
):
    assert taxa in [100, 500, 1000, 10000]
    if taxa == 10000:
        assert density == 0
    else:
        assert density in [20, 50, 75, 100]

    number_of_experiments = 30 if taxa != 10000 else 10

    for i in range(number_of_experiments):
        # source_file = f"bcd/data/SMIDGenOutgrouped/{taxa}/{density}/Source_Trees/ModelSourceTrees/smo.{i}.modelSourceTrees.tre"
        source_file = f"bcd/data/SMIDGenOutgrouped/{taxa}/{density}/Source_Trees/RaxML/smo.{i}.sourceTreesOG.tre"
        model_file = f"bcd/data/SMIDGenOutgrouped/{taxa}/{density}/Model_Trees/pruned/smo.{i}.modelTree.tre"
        if verbosity >= 1:
            print(f"Results for {i} ({source_file}):")
        run_methods(source_file, model_file, methods, verbosity=verbosity)


if __name__ == "__main__":
    # simulated_experiment(100, 20)
    # file = "data/superfine/500-taxa/20/sm_data.5" # Sup doesn't resolve for 5
    # file = "data/superfine/500-taxa/20/sm_data.10"

    # report(file + ".source_trees", file + ".model_tree", False)
    # taxa = 400
    # methods = [SCS, BCD]
    # for i in range(10):
    #     print(f"Results for {i}:")
    #     file = f"birth_death/{taxa}_taxa/{i}"
    #     # file = f"data/superfine/500-taxa/100/sm_data.{i}"
    #     run_methods(file + ".source_trees", file + ".model_tree", methods)

    # taxa = 500
    # density = 20
    methods = [SCS, BCD]
    # run_experiment_smidgen(taxa, density, methods)

    # run_experiment_super_triplets(75, 10, methods)

    taxa = 10000
    density = 0
    run_experiment_smidgen_og(taxa, density, methods, verbosity=2)

    # file = "data/superfine/100-taxa/20/sm_data.3"
    # report(file + ".source_trees", file + ".model_tree", False)

    # taxa = 1000
    # density = 100
    # for i in range(10):
    #     file = f"data/superfine/{taxa}-taxa/{density}/sm_data.{i}"
    #     print("DOING", file)
    #     report(file + ".source_trees", file + ".model_tree", False)

    # for i in range(10):  # range(10):
    #     print(f"Results for {i}:")
    #     file = f"simulated_trees/simulated_tree_data/100_taxa/{i}"
    #     # file = f"data/superfine/500-taxa/100/sm_data.{i}"
    #     report(file + ".source_trees", file + ".model_tree", False)

    # file = f"birth_death/100_taxa/9"
    # report(file + ".source_trees", file + ".model_tree", False)
    # scs = (
    #     make_tree(
    #         "(((((((((((t70,t46),(t55,(t39,t31),(t63,t86))),(((((t89,t26),t84,t48),((t42,t61,t33),(t99,(t34,t68,t62),t43))),(((t3,t51,t4),(t88,t50)),t75)),t1)),((t54,t81,t67),t96,t58)),(t77,t93)),((((((((t65,t14),t25,t28),(t76,t30,t9),t2),(((t45,t98),(t100,t87,t15),(t18,t95,t22)),((t19,t72),((t71,t79,t41),(t80,t69,t92))))),(t60,t5,t27)),((((t10,t85,t94),(t38,t12,t8)),(t56,t78,t16),t74),(t83,t11),t49)),(t36,(t13,t29),t90)),(t44,t64))),t23),(t91,t59)),(t20,(t21,t52,(t53,t40)),t73,(t17,t47))),t35),t97);"
    #     )
    #     .unrooted()
    #     .bifurcating()
    # ).unrooted()
    # sup = (
    #     make_tree(
    #         "(t94:0.031097998029898082,((t38:0.01587405493792574,(t12:0.1294589555996583,t8:0.02391769662423016):0.11893227254519337):0.027291540958539973,((t74:0.6387548273273306,(t56:0.6058782061160624,(t16:0.16878031781066216,t78:0.18930713301246826):0.4411361826834618):0.06632376577669119):0.13113858920592722,(((t36:0.5649249184402562,(t29:0.015284442584047313,t13:0.02584682185557445):0.05433861464524367):0.1798574394998838,(t90:0.6833691935583254,(t44:0.25762885506697525,t64:0.43127564558336856):0.3126971506190803)):0.17018434249289463,((t49:0.9235854549516707,(t83:0.3052683346288235,t11:0.46598028381208423):0.3326596129069099):0.31350967861611095,(((t73:0.12158970075746661,(t17:0.07812470280340002,t47:0.07342982946720784):0.03268263475270067):0.018555050565316228,(((t20:0.352193101930437,(t21:0.013808668586101558,(t53:0.05441672651732474,t40:0.007894214986670749):0.015237287529651145):0.05447582152619742):0.3217511518412046,(t35:0.4947064022532579,(t97:0.2252858055131954,t52:0.075113946570511):0.2261743014146011):0.30823816086623873):0.042429499201020515,((t91:0.16618820563780454,t59:0.31251893286372173):1.1656357497181717,(t23:1.1971748386498837,(t1:0.5998348603891359,((((((t61:0.02928022606290094,t42:0.016409349668629586):0.09158018077171343,t33:0.12850653967935433):0.21430834843080995,(t43:0.03997662031383534,((t62:0.05341894876476083,(t68:0.0672318306193603,t34:0.037877967953689785):0.03280626890716802):0.21923251669158927,t99:0.29737577826709394)):0.020022720430581165):0.03797618515640398,(((t26:0.2429774828010695,t89:0.07547669151171298):0.12429893547623495,t84:0.31540322403789184):0.01262282220290127,t48:0.4223212379437335):0.06338318923627118):0.02164248055656113,(((t88:0.0942163296551629,t50:0.22309860606112264):0.29637566459160497,((t51:0.113977398406318,t3:0.05928912652660677):0.3773899197410217,t4:0.3105290569997304)):0.06573637523524335,t75:0.5947984020752501)):0.009770139150162725,((t55:0.07043433918694744,((t46:0.046040935320862224,t70:0.03917175633965259):0.3634782842190702,((t86:0.04706979749549938,t63:0.04983500618633025):0.32730787163561836,(t31:0.17032928092312227,t39:0.14454285487056442):0.08503309977121408))):0.02840216495888445,((((t54:0.06867416503518521,t81:0.0972356236280755):0.028002902676154403,t67:0.08882145167183689):0.15408591226616056,t96:0.20915417315790605):0.1134182018009081,((t93:0.10834396031461196,t77:0.054860418058781286):0.1712448296181741,t58:0.08957824489661376):0.123675658739296):0.131735124031828):0.005047711466689991):0.005863810688114577):0.2201195137350506):0.2583566674101961):0.021329672822882714):0.009103440315449497):0.007689126899400675,((t27:0.2113722694571462,(t5:0.1628485309813183,t60:0.14595160220442646):0.06689401492253289):1.048690333114285,(((((t71:0.31758517816245574,(t79:0.3053006755725396,t41:0.06040520532254743)):0.14143823893140994,(t80:0.537967534761394,(t69:0.0016133782737070914,t92:1.12071732071908e-06):0.34437660446564594)):0.14076214250098332,(t72:0.04044532666199919,t19:0.8725364826036514):0.028294699731503475):0.09618192236909245,(((t98:0.09789937342832421,t45:0.08333364574861317):0.06051154449518597,(t15:0.07586541183237644,(t100:0.09443021517644011,t87:0.020588789894834442)):0.08465581698896266):0.05039573626858241,((t22:0.0354723549958756,t95:0.02722247255483302):0.024984783533903043,t18:0.04751744325346913):0.18962503841050166):0.538475569860528):0.005827403189976876,((t25:0.5965259622448638,(t28:0.610839939951756,(t14:0.059840250831387165,t65:0.13353450571964162):0.5852510942293071)):0.2628197290027999,(t2:0.9476493754708896,(t9:0.04916592519464276,(t76:0.08016995932744282,t30:0.06478556593549832):0.592410282420711):0.02661408188704391)):0.13646877292066467):0.012616373585455684):1.12071732071908e-06):0.035498950358952286):0.1670610093885456):0.02257473018098699):0.22693334556669878):0.02502419752558213,(t85:0.06519052295887567,t10:0.019372933546173473):0.25710849816141035);"
    #     )
    #     .unrooted()
    #     .bifurcating()
    # ).unrooted()
    # mcs = (
    #     make_tree(
    #         "((((t53,t20),((t40,t21),(t52,((((t47,t17,t73),((t64,(t44,(t90,(t72,((t60,t5,t27),((((t79,((t71,t80),(t69,t92))),t19),((t22,t18,t95),(t98,t45),(t87,t100))),((t76,(t28,(t14,t65),t25),t2),((t15,t41),(t30,((t11,t83,t49),((((t78,t16),t74,t56),((t94,t85,t10),(t38,t8,t12)),t9),(t29,t36,t13)))))))))))),(t1,(t77,(t93,((((t75,((t3,t51),(t50,t88,t4))),((t48,(t26,t89),t84),((t42,t33,t61),((t68,t62,t34),t99,t43)))),(t55,((t86,t63),((t31,t39),(t46,t70))))),((t54,t67,t81),t58,t96))))))),t23),(t91,t59))))),t35),t97);"
    #     )
    #     .unrooted()
    #     .bifurcating()
    # ).unrooted()
    # mod = (
    #     make_tree(
    #         "(((t17,t47),t73),((t91,t59),(((t23,((t20,(t21,(t53,t40))),((t97,t52),t35))),(t1,((t75,(((t3,t51),((t88,t50),t4)),((((t42,t61),t33),(t43,(((t68,t34),t62),t99))),((t84,(t89,t26)),t48)))),(((t46,t70),(((t86,t63),t55),(t39,t31))),(((t93,t77),t58),(t96,(t67,(t81,t54)))))))),(((((((t76,t30),t9),((t25,t28),(t14,t65))),t2),((((t22,t95),t18),((t100,(t15,t87)),(t45,t98))),((t19,t72),(((t79,t41),t71),((t69,t92),t80))))),((((t74,(t56,(t78,t16))),(((t8,t12),t38),(t94,(t85,t10)))),(t90,((t64,t44),((t13,t29),t36)))),(t49,(t11,t83)))),((t60,t5),t27)))));"
    #     )
    #     .unrooted()
    #     .bifurcating()
    # ).unrooted()

    # print(mod.lin_rajan_moret(scs))

    # print(mod.lin_rajan_moret(sup))

    # print(mod.lin_rajan_moret(mcs))
