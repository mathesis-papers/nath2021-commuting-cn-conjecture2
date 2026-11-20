import argparse
import json
from itertools import combinations, product
from pathlib import Path

from verify import characteristic_polynomial, common_matrix, polynomial_from_roots, polynomial_product


def dense_case(order, removed, parts, quotient, remaining_roots):
    adjacency = [[int(first != second) for second in range(order)] for first in range(order)]
    for first, second in removed:
        adjacency[first][second] = adjacency[second][first] = 0
    matrix = common_matrix(adjacency)
    assert sorted(vertex for part in parts for vertex in part) == list(range(order))
    for row, part in enumerate(parts):
        for vertex in part:
            assert [sum(matrix[vertex][neighbor] for neighbor in other) for other in parts] == quotient[row]
    quotient_polynomial = characteristic_polynomial(quotient)
    expected = polynomial_product(quotient_polynomial, polynomial_from_roots(remaining_roots))
    assert characteristic_polynomial(matrix) == expected
    assert len(remaining_roots) + len(parts) == order
    assert all(root < 0 for root in remaining_roots)
    if len(parts) == 2:
        assert quotient_polynomial[-1] < 0
    else:
        assert quotient_polynomial == [1, -19, -136, -144]
        assert -quotient_polynomial[-1] > 0 and quotient_polynomial[-2] < 0
    assert max(map(sum, matrix)) <= (order - 1) * (order - 2)
    return {'order': order, 'removed_edges': removed, 'part_sizes': list(map(len, parts)),
            'quotient': quotient, 'quotient_charpoly_descending': quotient_polynomial,
            'remaining_eigenvalues_with_multiplicity': remaining_roots,
            'full_charpoly_descending': expected, 'positive_eigenvalues': 1}


def dense_checks():
    output = []
    for order in [5, 6, 7]:
        quotient = [[order - 2, (order - 2) * (order - 3)],
                    [2 * (order - 3), (order - 2) * (order - 3)]]
        output.append(dense_case(order, [(0, 1)], [[0, 1], list(range(2, order))],
                                 quotient, [-(order - 2)] * (order - 2)))
    output.append(dense_case(7, [(0, 1), (2, 3)], [[0, 1, 2, 3], [4, 5, 6]],
                             [[11, 12], [16, 10]], [-5] * 4 + [-1]))
    output.append(dense_case(7, [(0, 1), (0, 2)], [[0], [1, 2], [3, 4, 5, 6]],
                             [[0, 8, 12], [4, 4, 16], [3, 8, 15]], [-4] + [-5] * 3))
    for order in range(3, 8):
        adjacency = [[int(first != second) for second in range(order)] for first in range(order)]
        matrix = common_matrix(adjacency)
        roots = [(order - 1) * (order - 2)] + [-(order - 2)] * (order - 1)
        assert characteristic_polynomial(matrix) == polynomial_from_roots(roots)
    thresholds = {order: ((2 * (order - 1) * (order - 2) + 2 * order) // 4)
                  for order in range(3, 8)}
    assert list(thresholds.values()) == [2, 5, 8, 13, 18]
    for order in range(8, 10001):
        assert 2 * (order - 1) * (order - 2) - 10 * order == 2 * (order - 4)**2 - 28
        assert 2 * (order - 4)**2 - 28 >= 4
    return {'sparse_edge_thresholds': thresholds,
            'complete_graph_orders': [3, 4, 5, 6, 7], 'noncomplete_dense_cases': output,
            'surface_polynomial_identity_integer_range': [8, 10000],
            'scope': 'Exact quotient matrices and characteristic polynomials support the small-order written proof; the all-order surface argument uses Euler and the written polynomial comparison.'}


def neighborhood_decompositions():
    count = 0
    for order in range(7):
        pairs = list(combinations(range(order), 2))
        for mask in range(1 << len(pairs)):
            adjacency = [[0] * order for _ in range(order)]
            for index, (first, second) in enumerate(pairs):
                if mask & (1 << index):
                    adjacency[first][second] = adjacency[second][first] = 1
            matrix = common_matrix(adjacency)
            decomposition = [[0] * order for _ in range(order)]
            for vertex in range(order):
                neighbors = [neighbor for neighbor in range(order) if adjacency[vertex][neighbor]]
                for first, second in combinations(neighbors, 2):
                    decomposition[first][second] += 1
                    decomposition[second][first] += 1
            assert matrix == decomposition
            degrees = list(map(sum, adjacency))
            assert 2 * sum(max(degree - 1, 0) for degree in degrees) == (
                2 * sum(degrees) - 2 * sum(degree > 0 for degree in degrees))
            count += 1
    return {'vertex_range': [0, 6], 'labelled_graphs': count,
            'scope': 'Entrywise neighborhood-clique decomposition and degree identity, not an exhaustive spectral proof.'}


def multiply(first, second, modulus=2):
    left_top, right_top, left_bottom, right_bottom = first
    other_left_top, other_right_top, other_left_bottom, other_right_bottom = second
    return ((left_top * other_left_top + right_top * other_left_bottom) % modulus,
            (left_top * other_right_top + right_top * other_right_bottom) % modulus,
            (left_bottom * other_left_top + right_bottom * other_left_bottom) % modulus,
            (left_bottom * other_right_top + right_bottom * other_right_bottom) % modulus)


def add(first, second, modulus=2):
    return tuple((left + right) % modulus for left, right in zip(first, second))


def ring_check(name, elements, modulus=2):
    positions = {element: index for index, element in enumerate(elements)}
    order = len(elements)
    addition = [[positions[add(first, second, modulus)] for second in elements] for first in elements]
    multiplication = [[positions[multiply(first, second, modulus)] for second in elements] for first in elements]
    assert elements[0] == (0, 0, 0, 0)
    for first, second, third in product(range(order), repeat=3):
        assert multiplication[multiplication[first][second]][third] == multiplication[first][multiplication[second][third]]
        assert multiplication[first][addition[second][third]] == addition[multiplication[first][second]][multiplication[first][third]]
        assert multiplication[addition[first][second]][third] == addition[multiplication[first][third]][multiplication[second][third]]
        assert addition[addition[first][second]][third] == addition[first][addition[second][third]]
    assert all(addition[first][0] == first and 0 in addition[first] for first in range(order))
    assert all(addition[first][second] == addition[second][first] for first, second in product(range(order), repeat=2))
    center = [first for first in range(order) if all(multiplication[first][second] == multiplication[second][first] for second in range(order))]
    vertices = [first for first in range(order) if first not in center]
    assert vertices
    vertex_count, center_order = len(vertices), len(center)
    adjacency = [[int(first != second and multiplication[first][second] == multiplication[second][first])
                  for second in vertices] for first in vertices]
    degrees = list(map(sum, adjacency))
    centralizers = [[second for second in range(order) if multiplication[first][second] == multiplication[second][first]] for first in vertices]
    for degree, centralizer in zip(degrees, centralizers):
        assert degree == len(centralizer) - center_order - 1
        assert order % len(centralizer) == 0 and 2 * len(centralizer) <= order
        assert all(addition[first][second] in centralizer for first, second in product(centralizer, repeat=2))
    remaining = set(range(vertex_count))
    components = []
    while remaining:
        pending = [min(remaining)]
        component = set(pending)
        while pending:
            current = pending.pop()
            for neighbor, adjacent in enumerate(adjacency[current]):
                if adjacent and neighbor not in component:
                    component.add(neighbor)
                    pending.append(neighbor)
        assert all(adjacency[first][second] == int(first != second) for first, second in product(component, repeat=2))
        remaining -= component
        components.append(len(component))
    roots = []
    for size in components:
        roots.extend([(size - 1) * (size - 2)] + [-(size - 2)] * (size - 1))
    matrix = common_matrix(adjacency)
    assert characteristic_polynomial(matrix) == polynomial_from_roots(roots)
    energy = sum(map(abs, roots))
    profile = 2 * sum(max(degree - 1, 0) for degree in degrees)
    bound = vertex_count * (vertex_count - center_order - 2)
    benchmark = 2 * (vertex_count - 1) * (vertex_count - 2)
    assert energy <= profile <= 2 * sum(degrees) <= bound < benchmark
    identities = [first for first in range(order) if all(multiplication[first][second] == second == multiplication[second][first] for second in range(order))]
    return {'name': name, 'order': order, 'has_multiplicative_identity': bool(identities),
            'center_order': center_order, 'noncentral_vertices': vertex_count,
            'centralizer_orders': list(map(len, centralizers)), 'degrees': degrees,
            'clique_component_sizes': sorted(components), 'charpoly_descending': characteristic_polynomial(matrix),
            'energy': energy, 'profile_bound': profile, 'centralizer_bound': bound, 'complete_graph_energy': benchmark}


def verify_extensions():
    matrices = list(product(range(2), repeat=4))
    rings = [ring_check('F2 row ring, matrices (a,b;0,0)', [entry for entry in matrices if entry[2] == entry[3] == 0]),
             ring_check('UT2(F2)', [entry for entry in matrices if entry[2] == 0]),
             ring_check('M2(F2)', matrices),
             ring_check('UT2(F3)', [entry for entry in product(range(3), repeat=4) if entry[2] == 0], modulus=3)]
    assert rings[0]['has_multiplicative_identity'] is False
    return {'schema_version': 1, 'status': 'pass', 'arithmetic': 'exact integers and rationals',
            'neighborhood_clique_decomposition': neighborhood_decompositions(),
            'small_order_and_surface_checks': dense_checks(), 'finite_ring_models': rings,
            'scope': 'Finite identities and models supplement the complete written proofs in note.tex.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    output = json.dumps(verify_extensions(), indent=2, ensure_ascii=False) + '\n'
    if arguments.output:
        arguments.output.write_text(output)
    else:
        print(output, end='')
