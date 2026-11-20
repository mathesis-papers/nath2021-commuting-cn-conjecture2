import argparse
from fractions import Fraction
from itertools import combinations, permutations
import json
from pathlib import Path


def multiply_matrices(left, right):
    columns = tuple(zip(*right))
    return [[sum(first * second for first, second in zip(row, column))
             for column in columns] for row in left]


def characteristic_polynomial(matrix):
    size = len(matrix)
    auxiliary = [[int(row == column) for column in range(size)]
                 for row in range(size)]
    coefficients = [1]
    for degree in range(1, size + 1):
        auxiliary = multiply_matrices(matrix, auxiliary)
        trace = sum(auxiliary[index][index] for index in range(size))
        assert trace % degree == 0
        coefficient = -trace // degree
        coefficients.append(coefficient)
        for index in range(size):
            auxiliary[index][index] += coefficient
    return coefficients


def polynomial_product(left, right):
    output = [0] * (len(left) + len(right) - 1)
    for first_index, first in enumerate(left):
        for second_index, second in enumerate(right):
            output[first_index + second_index] += first * second
    return output


def polynomial_from_roots(roots):
    result = [1]
    for root in roots:
        result = polynomial_product(result, [1, -root])
    return result


def common_matrix(adjacency):
    size = len(adjacency)
    literal = [[0 if row == column else sum(
        adjacency[row][vertex] * adjacency[column][vertex]
        for vertex in range(size)) for column in range(size)]
        for row in range(size)]
    square = multiply_matrices(adjacency, adjacency)
    for row in range(size):
        for column in range(size):
            assert literal[row][column] == square[row][column] - (
                sum(adjacency[row]) if row == column else 0)
    return literal


def all_graph_identities():
    count = 0
    for size in range(7):
        pairs = list(combinations(range(size), 2))
        for mask in range(1 << len(pairs)):
            adjacency = [[0] * size for _ in range(size)]
            for index, (first, second) in enumerate(pairs):
                if mask & (1 << index):
                    adjacency[first][second] = adjacency[second][first] = 1
            common_matrix(adjacency)
            count += 1
    return count


def partitions(total, minimum=1):
    if total == 0:
        yield ()
    for first in range(minimum, total + 1):
        for tail in partitions(total - first, first):
            yield (first,) + tail


def clique_checks():
    count = 0
    for size in range(1, 25):
        benchmark = 2 * (size - 1) * (size - 2)
        for components in partitions(size):
            energy = sum(2 * (part - 1) * (part - 2)
                         for part in components)
            assert energy <= benchmark
            assert (energy == benchmark) == (
                len(components) == 1 or components == (1, 1))
            count += 1
    return count


def symmetric_group(degree, alternating=False):
    elements = list(permutations(range(degree)))
    if alternating:
        elements = [entry for entry in elements if sum(
            entry[first] > entry[second]
            for first, second in combinations(range(degree), 2)) % 2 == 0]
    def multiplication(first, second):
        return tuple(first[index] for index in second)
    return elements, multiplication


def dihedral_group(rotation_order):
    elements = [(rotation, reflection) for rotation in range(rotation_order)
                for reflection in range(2)]
    def multiplication(first, second):
        return ((first[0] + (-1 if first[1] else 1) * second[0]) % rotation_order,
                (first[1] + second[1]) % 2)
    return elements, multiplication


def quaternion_group():
    elements = [(sign, basis) for sign in range(2) for basis in range(4)]
    basis_product = [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
    sign_product = [[0, 0, 0, 0], [0, 1, 0, 1], [0, 1, 1, 0], [0, 0, 1, 1]]
    def multiplication(first, second):
        return ((first[0] + second[0] + sign_product[first[1]][second[1]]) % 2,
                basis_product[first[1]][second[1]])
    return elements, multiplication


def group_check(name, elements, multiplication):
    positions = {element: index for index, element in enumerate(elements)}
    size = len(elements)
    table = [[positions[multiplication(first, second)] for second in elements]
             for first in elements]
    assert all(len(set(row)) == size for row in table)
    assert all(len(set(column)) == size for column in zip(*table))
    assert any(all(table[identity][index] == index == table[index][identity]
                   for index in range(size)) for identity in range(size))
    assert all(table[table[first][second]][third] == table[first][table[second][third]]
               for first in range(size) for second in range(size)
               for third in range(size))
    center = [first for first in range(size) if all(
        table[first][second] == table[second][first] for second in range(size))]
    vertices = [entry for entry in range(size) if entry not in center]
    order, center_order = len(vertices), len(center)
    assert order > 0
    adjacency = [[int(first != second and table[first][second] == table[second][first])
                  for second in vertices] for first in vertices]
    degrees = [sum(row) for row in adjacency]
    centralizer_orders = [sum(table[first][second] == table[second][first]
                              for second in range(size)) for first in vertices]
    assert all(degree + center_order + 1 == centralizer
               for degree, centralizer in zip(degrees, centralizer_orders))
    assert all(2 * centralizer <= size and size % centralizer == 0
               for centralizer in centralizer_orders)
    assert all(2 * degree <= order - center_order - 2 for degree in degrees)
    matrix = common_matrix(adjacency)
    polynomial = characteristic_polynomial(matrix)
    remaining = set(range(order))
    components = []
    while remaining:
        component = {min(remaining)}
        pending = list(component)
        while pending:
            current = pending.pop()
            for vertex, adjacent in enumerate(adjacency[current]):
                if adjacent and vertex not in component:
                    component.add(vertex)
                    pending.append(vertex)
        remaining -= component
        components.append(component)
    upper_bound = order * (order - center_order - 2)
    benchmark = 2 * (order - 1) * (order - 2)
    if name == "S4":
        expected = polynomial_product(polynomial_from_roots([0] * 8 + [3] * 2 + [-1] * 11),
                                      [1, -5, -30])
        assert polynomial == expected
        assert 25 < 145 < 13**2
        assert Fraction(17 + 13) <= 2 * sum(degrees) <= upper_bound < benchmark
        energy = "17 + sqrt(145)"
        component_sizes = None
    else:
        roots = []
        for component in components:
            assert all(adjacency[first][second] == int(first != second)
                       for first in component for second in component)
            component_order = len(component)
            roots.append((component_order - 1) * (component_order - 2))
            roots.extend([-(component_order - 2)] * (component_order - 1))
        assert polynomial == polynomial_from_roots(roots)
        energy = sum(abs(root) for root in roots)
        assert energy <= 2 * sum(degrees) <= upper_bound < benchmark
        component_sizes = sorted(map(len, components))
    return {"name": name, "group_order": size, "center_order": center_order,
            "vertex_count": order, "edge_count": sum(degrees) // 2,
            "common_neighborhood_charpoly_descending": polynomial,
            "clique_component_sizes": component_sizes, "exact_energy": energy,
            "centralizer_bound": upper_bound, "complete_graph_energy": benchmark}


def verify():
    group_models = [("S3", symmetric_group(3)), ("A4", symmetric_group(4, True)),
                    ("S4", symmetric_group(4)), ("Q8", quaternion_group())]
    group_models.extend((f"D{2 * rotation_order}", dihedral_group(rotation_order))
                        for rotation_order in range(3, 13))
    groups = [group_check(name, *model) for name, model in group_models]
    return {"schema_version": 1, "status": "pass", "arithmetic": "exact integers and rationals",
            "matrix_identity": {"vertex_range": [0, 6],
                                "labelled_graphs": all_graph_identities()},
            "clique_union_bound_and_equality": {"vertex_range": [1, 24],
                                               "integer_partitions": clique_checks()},
            "finite_group_models": groups,
            "scope": "Finite checks of the matrix identity, group centralizers, characteristic polynomials, and clique equality cases; universal conclusions have written proofs in note.tex."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    output = json.dumps(verify(), indent=2, ensure_ascii=False) + "\n"
    if arguments.output:
        arguments.output.write_text(output)
    else:
        print(output, end="")
