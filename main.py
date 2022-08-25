import fem

DEBUG = True


def main():
    node_pos = [(18.28, 9.14), (18.28, 0), (9.14, 9.14), (9.14, 0), (0, 9.14), (0, 0)]
    elem_node = [(2, 4), (0, 2), (3, 5), (1, 3), (2, 3), (0, 1), (3, 4), (2, 5), (1, 2), (0, 3)]

    sys = fem.FEM()
    for x, y in node_pos:
        sys.add_node(fem.Node(x, y))
    for s, e in elem_node[:6]:
        sys.add_element(fem.Elem(sys.nodes[s], sys.nodes[e], r=0.1))
    for s, e in elem_node[6:]:
        sys.add_element(fem.Elem(sys.nodes[s], sys.nodes[e], r=0.05))

    if DEBUG:
        print("L of element 1:", sys.elems[0].L)
        print("A of element 1:", sys.elems[0].A)
        print("K of element 1:", sys.elems[0].k)
        print("K of system:", sys.k / 1e6)


if __name__ == '__main__':
    main()
