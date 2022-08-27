from math import pi, sqrt, cos, sin
import cmath
import numpy as np
from typing import Union, Optional
from copy import copy

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

class Coord:
    """ Implementation of coordinate operation """

    def __init__(self, x, y):
        self._val = complex(x, y)

    @property
    def polar(self):
        return cmath.polar(self._val)

    @property
    def x(self):
        return self._val.real

    @property
    def y(self):
        return self._val.imag

    @property
    def ang(self):
        return self.polar[1]

    @property
    def cos(self):
        return cos(self.ang)

    @property
    def sin(self):
        return sin(self.ang)

    def __eq__(self, op):
        return self._val == op._val

    def __sub__(self, op):
        if not issubclass(type(op), Coord):
            raise TypeError("Type 'Coord' should only subtract type 'Coord'")
        return Coord(self.x - op.x, self.y - op.y)

    def __hash__(self):
        return hash(self._val)

    def __abs__(self):
        return abs(self._val)


class Node(Coord):
    ''' Class used for specifying nodes '''
    pass


class Elem:
    ''' Implementation of base class of elements used in FEM '''

    def __init__(self, node_i: Node, node_j: Node, r=0.1, E=2e11, sy=2.5e8, dens=7860):
        self.ps = node_i  # start position
        self.pe = node_j  # end position
        self._r = r  # radius
        self.E = E  # Young's module
        self.sy = sy  # yield stress
        self.dens = dens  #density

    r = property(lambda self: self._r)
    vec = property(lambda self: self.pe - self.ps)
    L = property(lambda self: abs(self.vec))
    A = property(lambda self: pi * self._r**2)
    ang = property(lambda self: self.vec.ang)
    cos = property(lambda self: self.vec.cos)
    sin = property(lambda self: self.vec.sin)

    @r.setter
    def r(self, val):
        self._r = val

    @property
    def trans_mat(self):
        ''' transfrom matrix from global coordinate to local one '''
        c = self.cos
        s = self.sin
        return np.array([c, s, -c, -s]).reshape((1, 4))

    @property
    def k(self):
        ''' stiffness matrix of the element '''
        c = self.cos
        s = self.sin
        tm = self.trans_mat
        k = (self.A * self.E / self.L) * np.dot(np.transpose(tm), tm)
        return k


class FEM:
    ''' class for FEM '''

    def __init__(self):
        self.nodes = []
        self.elems = []
        self.fix_nodes = []  # node indexes of fixed nodes
        self._F = None

    @property
    def F(self):
        ''' Forces applied on the truss '''
        return self._F

    @F.setter
    def F(self, F):
        assert len(F) == 2 * len(
            self.nodes), "Dimension of force applied doesn't match the number of nodes"
        self._F = F

    @property
    def fix_axis(self):
        ret = []
        for n in self.fix_nodes:
            ret += self.idx_n2a(n)
        return ret

    def node_idx(self, a: Union[float, Node], b: Optional[float] = None):
        '''
        Get the index of nodes in self.nodes. There are 2 ways to call this function
        
        `FEM.node_idx(node)` is used to search for a certain node
        `FEM.node_idx(x, y)` is used to search for a node with coordinate (x, y)
        '''

        tgt = a if b is None else Node(a, b)
        assert issubclass(type(tgt), Node), f"tgt is type {type(tgt)} but not Node"
        return self.nodes.index(tgt)

    def add_node(self, node: Node):
        ''' Add node to the FEM instance '''
        if node not in self.nodes:
            self.nodes.append(node)

    def add_element(self, elem: Elem):
        ''' Add element to the FEM instance '''
        if elem not in self.elems:
            self.elems.append(elem)

    def add_fix(self, node_idx):
        if node_idx not in self.fix_nodes:
            self.fix_nodes.append(node_idx)

    def idx_n2a(self, idx):
        ''' convert node index to axis index '''
        assert idx < len(self.nodes), f"{idx} is out of node indexes"
        return [2 * idx, 2 * idx + 1]

    def node_idx_of_elem(self, elem_idx):
        ''' Get the indexes of nodes (start & end) of the elements with given index '''
        idx = elem_idx
        return [self.node_idx(node) for node in [self.elems[idx].ps, self.elems[idx].pe]]

    def axis_idx_of_elem(self, elem_idx):
        ''' Get the indexes of axis of nodes (start & end) of the elements with given index '''
        idx = elem_idx
        n_idxs = self.node_idx_of_elem(idx)
        return self.idx_n2a(n_idxs[0]) + self.idx_n2a(n_idxs[1])

    @property
    def k(self):
        ''' Stiffness matrix of the whole system '''
        size = 2 * len(self.nodes)
        k = np.zeros((size, size))

        for i in range(len(self.elems)):
            elem = self.elems[i]
            axis_idxs = self.axis_idx_of_elem(i)
            l = len(axis_idxs)
            for r in range(l):
                for c in range(l):
                    k[axis_idxs[r]][axis_idxs[c]] += elem.k[r][c]
        return k

    @property
    def axis_displacement(self):
        ''' displacement of every axis when force applied '''
        assert self.F, "Force are not applied"
        size = 2 * len(self.nodes)
        k = copy(self.k)
        F = copy(self.F)
        idxes = np.array(range(size))

        idxes = np.delete(idxes, self.fix_axis, None)
        k = np.delete(k, self.fix_axis, 0)
        k = np.delete(k, self.fix_axis, 1)
        F = np.delete(F, self.fix_axis, 0)

        k_inv = np.linalg.inv(k)
        Q = np.dot(k_inv,F)
        out = np.zeros(size)
        for i, idx in zip(range(len(idxes)), idxes):
            out[idx] = Q[i]
        return out

    @property
    def node_displacement(self):
        a_disp = self.axis_displacement
        disp = np.zeros(len(self.nodes))
        for i , val in zip(range(len(disp)), chunks(a_disp, 2)):
            disp[i] = abs(complex(*val))
        return disp

    @property
    def stress(self):
        ''' stress of each elements when force applied '''
        stress = []
        for i in range(len(self.elems)):
            e = self.elems[i]
            axis_idxs = self.axis_idx_of_elem(i)
            disp = np.array([self.axis_displacement[i] for i in axis_idxs])
            s = e.E / e.L * np.dot(e.trans_mat, disp)
            stress.append(s)
        return stress

    @property
    def total_mass(self):
        return sum([e.dens * e.A * e.L for e in self.elems])

    def over_loading(self):
        ''' difference between applied stress and yielding stress for all elements '''
        strs = self.stress
        return np.array([abs(strs[i]) - self.elems[i].sy for i in range(len(self.elems))])

    def set_elem_rad(self, elem_idx, radius):
        self.elems[elem_idx].r = radius
