from quiver import *
# from moduli import *  # ?

quiver_data_light = [
    [GeneralizedKroneckerQuiver(5), vector([2, 3]), vector([3, -2])],
    # [SubspaceQuiver(5), [1 for i in range(5)] + [5], [-1 for i in range(5)] + [1]],
    # [ThreeVertexQuiver(1, 2, 1), [2, 1, 3], [1, 1, -1]]
]

quiver_data_heavy = [
    [GeneralizedKroneckerQuiver(17), vector([7, 13]), vector([13, -7])],
    # [SubspaceQuiver(10), [1 for i in range(10)] + [10], [-1 for i in range(10)] + [1]],
    # [ThreeVertexQuiver(1, 6, 7), [4, 1, 4], [25, 24, -31]]
]

batch_0_data = [
    [GeneralizedKroneckerQuiver(5)],
    # [SubspaceQuiver(20)],
    # [ThreeVertexQuiver(1, 6, 7)],
]

results = {}

# testing is_acyclic
results["is_acyclic()"] = {}
results["is_connected()"] = {}

for datum in batch_0_data:
    datum[0].is_acyclic()
    results["is_acyclic()"][datum[0]._name] = timeit("datum[0].is_acyclic()", number=1000)
    results["is_connected()"][datum[0]._name] = timeit("datum[0].is_connected()", number=1000)

# testing all_HN_types
results["all_harder_narasimhan_types()"] = {}
# results["all_weight_bounds()"] = {}
results["euler_form()"] = {}
results["has_semistable_representation()"] = {}

for datum in quiver_data_light + quiver_data_heavy:
    Q = datum[0]
    d = datum[1]
    theta = datum[2]
    M = QuiverModuliStack(Q, d, theta)
    results["all_harder_narasimhan_types()"][Q._name] = timeit(
        "M.all_harder_narasimhan_types()", number=100
    )
    # results["all_weight_bounds()"][Q._name] = timeit(
    #     "Q.all_weight_bounds(d, theta)", number=100
    # )
    results["euler_form()"][Q._name] = timeit("Q.euler_form(d, d)", number=1000)
    results["has_semistable_representation()"][Q._name] = timeit(
        "Q.has_semistable_representation(d, theta)", number=1000
    )

# results["generic_ext()"] = {}
# results["generic_hom()"] = {}
results["is_schur_root()"] = {}
results["canonical_decomposition()"] = {}
results["in_fundamental_domain()"] = {}

for datum in quiver_data_light + quiver_data_heavy:
    Q = datum[0]
    d = datum[1]
    theta = datum[2]
    M = QuiverModuli(datum[0], datum[1], datum[2])
    # results["generic_ext()"][datum[0]] = timeit(
    #     "M.generic_ext(d, d)", number=100
    # )
    # results["generic_hom()"][datum[0]] = timeit(
    #     "M.generic_hom(d, d)", number=100
    # )
    results["is_schur_root()"][datum[0]] = timeit(
        "M.is_schur_root(d)", number=100
    )
    results["canonical_decomposition()"][datum[0]] = timeit(
        "M.canonical_decomposition(d)", number=100
    )
    results["in_fundamental_domain()"][datum[0]] = timeit(
        "M.in_fundamental_domain(d)", number=100
    )


# save results to file
with open("benchmarks_results.md", "w") as f:
    for key in results:
        f.write(f"## {key}\n")
        f.write("| Quiver | Time |\n")
        f.write("| --- | --- |\n")
        for q in results[key]:
            f.write(f"| {q} | {results[key][q]} |\n")
        f.write("\n")
        f.write("\n")
