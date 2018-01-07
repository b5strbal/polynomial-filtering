from itertools import chain, combinations
 
def all_subsets(ss):
    """
    Return the list of all subsets of a set.
    """
    return chain(*map(lambda x: combinations(ss, x), range(0, len(ss)+1)))
 
# for subset in all_subsets([1, 2, 3, 4]):
    #   print(subset)

class SymmetricGraphPA(object):
    """
    A pseudo-Anosov mapping class constructed using Penner's construction for a rotationally symmetric graph as a composition of a twist and rotation.
    """
    def __init__(self, num_verts, connected_distances,
                outgoing_edge_permutation=None):
        """
        Initialize the object.

        INPUT:

        - ``num_verts`` -- the number of vertices of the graph (equivalently, the number of curves in Penner's construction), arranged along the vertices of regular ``num_verts``-gon.
        - `connected_distances`` -- a list or set of positive integers, specifying the diagonals of the regular polygon that are part of the graph. I.e. if this set is {1,3}, then even side of the polygon is part of the graph and vertices of distance three are also connected.

        """
        if num_verts < 3:
            raise ValueError("The number of vertices has to be at least 3")
        graph = Graph(num_verts) # vertices numbered from 0 to num_verts
        for distance in connected_distances:
            if distance < 1 or distance >= num_verts:
                raise ValueError("The distances have to be between 1 and half the number of vertices minus 1 .")
            for v in graph.vertices():
                graph.add_edge(v, (v+distance) % num_verts)

        if not graph.is_connected():
            raise ValueError("The graph is not connected!")

        if outgoing_edge_permutation is None:
            # Default is the trivial permutation
            self._outgoing_edge_permutation = range(graph.degree(0))
        else:
            if sorted(outgoing_edge_permutation) != range(graph.degree(0)):
                raise ValueError("Invalid permutation.")
            else:
                self._outgoing_edge_permutation = outgoing_edge_permutation
        self._graph = graph
        self._num_verts = num_verts
        self._connected_distances = set(connected_distances)

    def __repr__(self):
        return "SymmetricGraphPA({0},{1})".format(self._num_verts, list(self._connected_distances))

    @classmethod
    def with_consecutive_neighbors(cls, num_verts, valence):
        """
        Construct an example where the neighbors of each vertex are consecutive.

        INPUT:

        - ``num_verts`` -- the number of vertices
        - ``valence`` -- the valence of each vertex

        EXAMPLE:

        sage: pa = SymmetricGraphPA.with_consecutive_neighbors(7,2)
        sage: pa.poly()
        x^7 - x^4 - x^3 - 1

        sage: pa = sage: SymmetricGraphPA.with_consecutive_neighbors(8,3)
        sage: pa.poly()
        x^8 - x^5 - x^4 - x^3 - 1

        """
        if num_verts % 2 == valence % 2:
            raise ValueError("The parity of the number of vertices and valence have to be different.")
        return cls(num_verts, range((num_verts+1)//2-valence//2, num_verts//2+1))

    def is_surface_orientable(self):
        """
        Decide if the surface is orientable.

        EXAMPLES:

            sage: SymmetricGraphPA(3, {2}).is_surface_orientable()
            False
            sage: SymmetricGraphPA(4, {1}).is_surface_orientable()
            True
        """
        return self._graph.is_bipartite()

    def poly(self):
        """
        Return a polynomial whose largest root is the stretch factor.

        The polynomials is the characteristic polynomial of the product of Penner matrices. This matrix might be bigger (maybe also smaller) than the action on homology.

        EXAMPLES:

            sage: SymmetricGraphPA(3, {1}).poly()
            x^3 - x^2 - x - 1
            sage: SymmetricGraphPA(7, {1,3}).poly()
            x^7 - x^6 - x^4 - x^3 - x - 1

        """
        num_verts = self._num_verts
        pol = x^num_verts - 1
        symmetric_distances = self._connected_distances.union(
            {num_verts-dist for dist in self._connected_distances})
        for i in symmetric_distances:
            pol -= x^i
        return pol

    def stretch_factor(self):
        """
        Return the stretch factor of the pseudo-Anosov map.

        EXAMPLE:

            sage: SymmetricGraphPA(3, {1}).stretch_factor()
            1.839286755214161?

        """
        roots = self.poly().roots(ring=QQbar)
        return max(abs(item[0]) for item in roots)

    def _next_neighbor(self, start_vertex, end_vertex, reverse=False):
        """
        Return the next neighbor of vertex in the Penner graph.

        This makes sense, because the there is a natural cyclic orientation of the neighbors inherited from their position along the regular polygon.

        INPUT:

        - ``start_vertex`` -- the vertex whose neighbors are considered
        - ``end_vertex`` -- a neighbor of ``start_vertex``
        - ``reverse`` -- (default:False) if True, the previous neighbor 
        is returned, not the next one.

        EXAMPLE:

            sage: pa = SymmetricGraphPA(6, {2,3})
            sage: pa._next_neighbor(1, 3)
            5
            sage: pa._next_neighbor(1, 4)
            3
            sage: pa._next_neighbor(1, 5)
            4
            sage: pa._next_neighbor(1, 3, True)
            4

        """
        while True:
            end_vertex += 1 if reverse else -1
            end_vertex %= self._num_verts
            if self._graph.has_edge(start_vertex, end_vertex):
                return end_vertex

    def _bdy_graph(self):
        """
        Return a graph consisting of cycles whose connected components are the boundary components of the surface.
        """
        bdy_graph = Graph()
        for start_vertex in self._graph.vertices():
            # consider the annulus corresponding to start_vertex
            for end_vertex in self._graph.neighbors(start_vertex):
                # consider the intersection of this annulus with the annulus
                # for end_vertex and also the intersection with the next
                # annulus
                for side in ['l', 'r']:
                    # finally, consider the interval between the two
                    # intersections on the left or right
                    # 1. forward neighbor
                    new_start_vertex = self._next_neighbor(start_vertex,
                                                            end_vertex)
                    if side == 'l':
                        new_end_vertex = self._next_neighbor(new_start_vertex,
                                                            start_vertex, True)
                    else:
                        new_end_vertex = start_vertex
                    bdy_graph.add_edge((start_vertex,end_vertex,side),
                                 (new_start_vertex, new_end_vertex, 'l'))
                    # 2. backward neighbor
                    if side == 'l':
                        new_end_vertex = self._next_neighbor(
                            end_vertex, start_vertex, True)
                    else:
                        new_end_vertex = start_vertex

                    bdy_graph.add_edge((start_vertex,end_vertex,side),
                                 (end_vertex, new_end_vertex, 'r'))
        return bdy_graph

    def euler_char(self):
        """
        Return the Euler characteristic of the surface with boundries.

        Bigon complementary regions are also considered boundaries.

        EXAMPLE:

            sage: SymmetricGraphPA(3, {1}).euler_char()
            -3
            sage: SymmetricGraphPA(22, {11,10,9,8,7,6}).euler_char()
            -121

        """
        return -self._graph.num_edges()

    def singularity_type(self):
        """
        Return the singularity type of the pseudo-Anosov.

        EXAMPLE:

            sage: SymmetricGraphPA(3, {1}).singularity_type()
            [6]
            sage: SymmetricGraphPA(22, {11,10,9,8,7,6}).singularity_type()
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

        """
        bdy_graph = self._bdy_graph()
        result = []
        for x in bdy_graph.connected_components():
            if len(x) > 4: # cycles of length 4 are removable singularities
                result.append(len(x)/2)
        return result       

    def num_boundaries(self):
        """
        Return the number of boundaries, including bigons.

        EXAMPLE:

            sage: SymmetricGraphPA(3, {1}).num_boundaries()
            1
            sage: SymmetricGraphPA(22, {11,10,9,8,7,6}).num_boundaries()
            110

        """
        # valence = self._graph.degree(0)
        bdy_graph = self._bdy_graph()
        count = 0
        for x in bdy_graph.connected_components():
            # if len(x) > 4: # cycles of length 4 are removable singularities
            count += 1
        return count
        # return gcd(self._num_verts, self._valence)

    def genus(self):
        """
        Return the genus of the surface.

        EXAMPLE:

            sage: SymmetricGraphPA(3, {1}).genus()
            4
            sage: SymmetricGraphPA(22, {11,10,9,8,7,6}).genus()
            13

        """
        if self.is_surface_orientable():
            # euler_char = 2-2*genus-num_boundaries
            twice_genus = 2-self.euler_char()-self.num_boundaries()
            assert(twice_genus % 2 == 0)
            return twice_genus // 2
        else:
            # euler_char = 2-genus-num_boundaries
            return 2-self.euler_char()-self.num_boundaries()


def generate_examples(num_verts_list, only_consecutive=True):
    """
    Compute the surfaces and stretch factors for all graphs with number of vertices in the specified set.

    INPUT:
    - ``num_verts_list`` -- the list of values that are tried as number of vertices
    - ``only_consecutive`` -- (default:True) if True, the graph are constructed only using the ``with_consecutive_neighbors`` constructor. This is fast and gives the best results for nonorientable surfaces. If False, all possible graph are tried. This get exponentially slow in the number of vertices, but constructs more examples, including orientation-reversing examples for orientable surfaces.

    EXAMPLES: 

        sage: generate_examples(range(5))
        {(False, 4): [(3, {1}, x^3 - x^2 - x - 1, 1.839286755214161?)],
        (False, 5): [(4, {1, 2}, x^4 - x^3 - x^2 - x - 1, 1.927561975482926?)],
        (False, 6): [(5, {2}, x^5 - x^3 - x^2 - 1, 1.429108319838146?),
        (5, {1, 2}, x^5 - x^4 - x^3 - x^2 - x - 1, 1.965948236645486?)]}

    To list the smallest stretch factors found for each surface, run

        sage: examples = generate_examples(range(10))
        sage: {key: examples[key][0] for key in examples}

        {(False, 4): (3, {1}, x^3 - x^2 - x - 1, 1.839286755214161?),
        (False, 5): (6, {2, 3}, x^6 - x^4 - x^3 - x^2 - 1, 1.512876396864095?),
        (False, 6): (5, {2}, x^5 - x^3 - x^2 - 1, 1.429108319838146?),
        (False, 7): (6,
        {1, 2, 3},
        x^6 - x^5 - x^4 - x^3 - x^2 - x - 1,
        1.983582843424327?),
        (False, 8): (7, {3}, x^7 - x^4 - x^3 - 1, 1.288452726276414?),
        (False, 9): (8, {3, 4}, x^8 - x^5 - x^4 - x^3 - 1, 1.356797155088499?),
        (False, 10): (9, {4}, x^9 - x^5 - x^4 - 1, 1.217281181744868?)} 

    """

    def try_example(pa):
        try:
            # print num_verts, connected_distances
            key = (pa.is_surface_orientable(), pa.genus())
            if key not in result:
                result[key] = []
            result[key].append((pa._num_verts, pa._connected_distances, pa.poly(), pa.stretch_factor()))
        except ValueError:
            pass
    result = {}
    for num_verts in num_verts_list:
        if only_consecutive:
            for valence in range(1, num_verts):
                try:
                    pa = SymmetricGraphPA.with_consecutive_neighbors(num_verts, valence)
                    try_example(pa)
                except ValueError:
                    pass
        else:
            distances = range(1,num_verts//2+1)
            for connected_distances in all_subsets(distances):
                pa = SymmetricGraphPA(num_verts, connected_distances)
                try_example(pa)
    for key in result:
        result[key].sort(key=lambda x: x[3])
    return result

# def find_examples_with_genus(g):
#     ret = []
#     for num_verts in range(3, 2*(g-2)+1):
#         for valence in range(num_verts-1,1,-2):
#             pa = SymmetricGraphPA(num_verts, valence)
#             if pa.genus() == g:
#                 ret.append((num_verts,valence,pa.stretch_factor()))
#     return ret

# def minimal_stretch_factors(min_g, max_g):
#     """
#     EXAMPLE:

#     sage: minimal_stretch_factors(4, 20)

#     """
#     ret = {}
#     for g in range(min_g,max_g+1):
#         ret[g] = min(find_examples_with_genus(g), key=lambda x: x[2])
#     return ret
