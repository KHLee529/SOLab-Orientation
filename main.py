import fem
from scipy import optimize
import numpy as np


class fem_opt:

    def __init__(self, sys):
        assert len(sys.elems) == 10
        self.sys = sys

    def set_rad(self, r1, r2):
        for i in range(6):
            self.sys.set_elem_rad(i, r1)
        for i in range(6, 10):
            self.sys.set_elem_rad(i, r2)

    def obj_func(self, x, *args):
        ''' objective function for optimization '''
        assert len(x) == 2, "x should be a tuple (r1, r2) with length 2"
        self.set_rad(*x)
        return self.sys.total_mass

    def disp_cons(self, x):
        self.set_rad(*x)
        return self.sys.node_displacement[1]

    def stress_cons(self, x):
        self.set_rad(*x)
        return self.sys.over_loading()


def main():
    node_pos = [(18.28, 9.14), (18.28, 0), (9.14, 9.14), (9.14, 0), (0, 9.14), (0, 0)]
    elem_node = [(2, 4), (0, 2), (3, 5), (1, 3), (2, 3), (0, 1), (3, 4), (2, 5), (1, 2), (0, 3)]
    F = [0, 0, 0, 1e7, 0, 0, 0, 1e7, 0, 0, 0, 0]

    sys = fem.FEM()
    for x, y in node_pos:
        sys.add_node(fem.Node(x, y))
    for s, e in elem_node[:6]:
        sys.add_element(fem.Elem(sys.nodes[s], sys.nodes[e], r=0.5))
    for s, e in elem_node[6:]:
        sys.add_element(fem.Elem(sys.nodes[s], sys.nodes[e], r=0.5))
    for i in range(4, 6):
        sys.add_fix(i)

    sys.F = F
    opt = fem_opt(sys)

    x0 = (0.5, 0.5)
    bnds = ((0.001, 0.5), (0.001, 0.5))
    cons = [
        optimize.NonlinearConstraint(opt.disp_cons, lb=-np.inf, ub=0.02),
    ] + [
        optimize.NonlinearConstraint(lambda x: opt.stress_cons(x)[i], lb=-np.inf, ub=0)
        for i in range(len(sys.elems))
    ]

    res = optimize.minimize(opt.obj_func, x0, method='SLSQP', bounds=bnds, constraints=cons)
    print(res)


if __name__ == '__main__':
    main()
