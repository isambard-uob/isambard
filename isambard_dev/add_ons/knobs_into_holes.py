import numpy
import networkx
import itertools
from collections import OrderedDict, Counter
from scipy.cluster.hierarchy import linkage, fcluster

from ampal.protein import Polypeptide
from ampal.assembly import Assembly
from ampal.pdb_parser import convert_pdb_to_ampal
from ampal.pseudo_atoms import PseudoAtom, PseudoGroup, PseudoMonomer, Primitive
from ampal.ampal_databases import element_data
from ampal.analyse_protein import polypeptide_vector, crick_angles
from .pacc import fit_heptad_register
from tools.geometry import centre_of_mass, distance, angle_between_vectors, is_acute, find_foot, Axis, \
    minimal_distance_between_lines



_heptad_colours = {
    'a': 'red',
    'b': 'orange',
    'c': 'yellow',
    'd': 'green',
    'e': 'cyan',
    'f': 'blue',
    'g': 'violet'
}


def side_chain_centres(assembly, masses=False):
    """ PseudoGroup containing side_chain centres of each Residue in each Polypeptide in Assembly.

    Notes
    -----
    Each PseudoAtom is a side-chain centre.
    There is one PseudoMonomer per chain in ampal (each containing len(chain) PseudoAtoms).
    The PseudoGroup has len(ampal) PseudoMonomers.

    Parameters
    ----------
    assembly : Assembly
    masses : bool
        If True, side-chain centres are centres of mass.
        If False, side-chain centres are centres of coordinates.

    Returns
    -------
    PseudoGroup
        containing all side_chain centres, and with ampal_parent=assembly.
    """
    if masses:
        elts = set([x.element for x in assembly.get_atoms()])
        masses_dict = {e: element_data[e]['atomic mass'] for e in elts}
    pseudo_monomers = []
    for chain in assembly:
        if isinstance(chain, Polypeptide):
            centres = OrderedDict()
            for r in chain.get_monomers(ligands=False):
                side_chain = r.side_chain
                if masses:
                    masses_list = [masses_dict[x.element] for x in side_chain]
                else:
                    masses_list = None
                if side_chain:
                    centre = centre_of_mass(points=[x._vector for x in side_chain], masses=masses_list)
                # for Glycine residues.
                else:
                    centre = r['CA']._vector
                centres[r.unique_id] = PseudoAtom(coordinates=centre, name=r.unique_id, ampal_parent=r)
            pseudo_monomers.append(PseudoMonomer(pseudo_atoms=centres, monomer_id=' ', ampal_parent=chain))
    return PseudoGroup(monomers=pseudo_monomers, ampal_parent=assembly)


def cluster_helices(helices, cluster_distance=12.0):
    """ Clusters helices according to the minimum distance between the line segments representing their backbone.

    Notes
    -----
    Each helix is represented as a line segement joining the CA of its first Residue to the CA if its final Residue.
    The minimal distance between pairwise line segments is calculated and stored in a condensed_distance_matrix.
    This is clustered using the 'single' linkage metric
     (all members of cluster i are at < cluster_distance away from at least one other member of cluster i).
    Helices belonging to the same cluster are grouped together as values of the returned cluster_dict.

    Parameters
    ----------
    helices: Assembly
    cluster_distance: float

    Returns
    -------
    cluster_dict: dict
        Keys: int
            cluster number
        Values: [Polymer]

    """
    condensed_distance_matrix = []
    for h1, h2 in itertools.combinations(helices, 2):
        md = minimal_distance_between_lines(h1[0]['CA']._vector, h1[-1]['CA']._vector,
                                            h2[0]['CA']._vector, h2[-1]['CA']._vector, segments=True)
        condensed_distance_matrix.append(md)
    z = linkage(condensed_distance_matrix, method='single')
    clusters = fcluster(z, t=cluster_distance, criterion='distance')
    cluster_dict = {}
    for h, k in zip(helices, clusters):
        if k not in cluster_dict:
            cluster_dict[k] = [h]
        else:
            cluster_dict[k].append(h)
    return cluster_dict


def find_kihs(assembly, hole_size=4, cutoff=7.0):
    """ KnobIntoHoles between residues of different chains in assembly.

    Notes
    -----
    A KnobIntoHole is a found when the side-chain centre of a Residue a chain is close than (cutoff) Angstroms from at
    least (hole_size) side-chain centres of Residues of a different chain.

    Parameters
    ----------
    assembly : Assembly
    hole_size : int
        Number of Residues required to form each hole.
    cutoff : float
        Maximum distance between the knob and each of the hole residues.

    Returns
    -------
    kihs : [KnobIntoHole]
    """
    pseudo_group = side_chain_centres(assembly=assembly, masses=False)
    pairs = itertools.permutations(pseudo_group, 2)
    kihs = []
    for pp_1, pp_2 in pairs:
        for r in pp_1:
            close_atoms = pp_2.is_within(cutoff, r)
            # kihs occur between residue and (hole_size) closest side-chains on adjacent polypeptide.
            if len(close_atoms) < hole_size:
                continue
            elif len(close_atoms) > hole_size:
                close_atoms = sorted(close_atoms, key=lambda x: distance(x, r))[:hole_size]
            kih = OrderedDict()
            kih['k'] = r
            for i, hole_atom in enumerate(close_atoms):
                kih['h{0}'.format(i)] = hole_atom
            knob_into_hole = KnobIntoHole(pseudo_atoms=kih)
            kihs.append(knob_into_hole)
    return kihs


def find_contiguous_packing_segments(polypeptide, residues, max_dist=10.0):
    """ Assembly containing segments of polypeptide, divided according to separation of contiguous residues.

    Parameters
    ----------
    polypeptide : Polypeptide
    residues : iterable containing Residues
    max_dist : float
        Separation beyond which splitting of Polymer occurs.

    Returns
    -------
    segments : Assembly
        Each segment contains a subset of residues, each not separated by more than max_dist from the previous Residue.
    """
    segments = Assembly(assembly_id=polypeptide.ampal_parent.id)
    residues_in_polypeptide = list(sorted(residues.intersection(set(polypeptide.get_monomers())),
                                          key=lambda x: int(x.id)))
    if not residues_in_polypeptide:
        return segments
    # residue_pots contains separate pots of residues divided according to their separation distance.
    residue_pots = []
    pot = [residues_in_polypeptide[0]]
    for r1, r2 in zip(residues_in_polypeptide, residues_in_polypeptide[1:]):
        d = distance(r1['CA'], r2['CA'])
        if d <= max_dist:
            pot.append(r2)
            if sum([len(x) for x in residue_pots] + [len(pot)]) == len(residues_in_polypeptide):
                residue_pots.append(pot)
        else:
            residue_pots.append(pot)
            pot = [r2]
    for pot in residue_pots:
        segment = polypeptide.get_slice_from_res_id(pot[0].id, pot[-1].id)
        segment.ampal_parent = polypeptide.ampal_parent
        segments.append(segment)
    return segments


class KnobGroup(PseudoGroup):
    """ Collection of KnobsIntoHoles interactions with associated methods. """
    def __init__(self, monomers=None, polymer_id=' ', ampal_parent=None, cutoff=None):
        super(KnobGroup, self).__init__(monomers=monomers, polymer_id=polymer_id, ampal_parent=ampal_parent)
        if cutoff is None:
            cutoff = max([x.max_kh_distance for x in self])
        self.cutoff = cutoff

    def __repr__(self):
        return '<KnobGroup containing {} {}>'.format(
            len(self._monomers), 'KnobIntoHole' if len(self._monomers) == 1 else 'KnobsIntoHoles')

    @classmethod
    def from_helices(cls, assembly, cutoff=7.0, min_helix_length=8):
        """ Generate KnobGroup from the helices in the assembly - classic socket functionality.

        Notes
        -----
        Socket identifies knobs-into-holes (KIHs) packing motifs in protein structures.
        The following resources can provide more information:
        The socket webserver: http://coiledcoils.chm.bris.ac.uk/socket/server.html
        The help page: http://coiledcoils.chm.bris.ac.uk/socket/help.html
        The original publication reference: Walshaw, J. & Woolfson, D.N. (2001) J. Mol. Biol., 307 (5), 1427-1450.

        Parameters
        ----------
        assembly : Assembly
        cutoff : float
            Socket cutoff in Angstroms
        min_helix_length : int
            Minimum number of Residues in a helix considered for KIH packing.

        Returns
        -------
        instance : KnobGroup
        None if no helices or no kihs.
        """
        cutoff = float(cutoff)
        helices = Assembly([x for x in assembly.helices if len(x) >= min_helix_length])
        if len(helices) <= 1:
            return None
        # reassign ampal_parents
        helices.relabel_polymers([x.ampal_parent.id for x in helices])
        for i, h in enumerate(helices):
            h.number = i
            h.ampal_parent = h[0].ampal_parent
            for r in h.get_monomers():
                r.tags['helix'] = h
        all_kihs = []
        cluster_dict = cluster_helices(helices, cluster_distance=(cutoff + 10))
        for k, v in cluster_dict.items():
            if len(v) > 1:
                kihs = find_kihs(v, cutoff=cutoff, hole_size=4)
                if len(kihs) == 0:
                    continue
                for x in kihs:
                    all_kihs.append(x)
        instance = cls(ampal_parent=helices, cutoff=cutoff)
        for x in all_kihs:
            x.ampal_parent = instance
        instance._monomers = all_kihs
        instance.relabel_monomers()
        return instance

    def knob_subgroup(self, cutoff=7.0):
        """ KnobGroup where all KnobsIntoHoles have max_kh_distance <= cutoff. """
        if cutoff > self.cutoff:
            raise ValueError("cutoff supplied ({0}) cannot be greater than self.cutoff ({1})".format(cutoff,
                                                                                                     self.cutoff))
        return KnobGroup(monomers=[x for x in self.get_monomers()
                                   if x.max_kh_distance <= cutoff], ampal_parent=self.ampal_parent)

    @property
    def complementary_knobs(self, cutoff=None):
        return list(set([x.knob_residue for x in self.get_monomers() if x.is_complementary(cutoff=cutoff)]))

    @property
    def graph(self):
        """ Returns MultiDiGraph from kihs. Nodes are helices and edges are kihs. """
        g = networkx.MultiDiGraph()
        edge_list = [(x.knob_helix, x.hole_helix, x.id, {'kih': x}) for x in self.get_monomers()]
        g.add_edges_from(edge_list)
        return g

    @staticmethod
    def filter_graph(g, cutoff=7.0, min_kihs=2):
        """ Get subgraph formed from edges that have max_kh_distance < cutoff.

        Parameters
        ----------
        g : MultiDiGraph representing KIHs
            g is the output from graph_from_protein
        cutoff : float
            Socket cutoff in Angstroms.
            Default is 7.0.
        min_kihs : int
            Minimum number of KIHs shared between all pairs of connected nodes in the graph.

        Returns
        -------
        networkx.MultiDigraph
            subgraph formed from edges that have max_kh_distance < cutoff.
        """
        edge_list = [e for e in g.edges(keys=True, data=True) if e[3]['kih'].max_kh_distance <= cutoff]
        if min_kihs > 0:
            c = Counter([(e[0], e[1]) for e in edge_list])
            # list of nodes that share > min_kihs edges with at least one other node.
            node_list = set(list(itertools.chain.from_iterable([k for k, v in c.items() if v > min_kihs])))
            edge_list = [e for e in edge_list if (e[0] in node_list) and (e[1] in node_list)]
        return networkx.MultiDiGraph(edge_list)

    def get_coiledcoil_region(self, cc_number=0, cutoff=7.0, min_kihs=2):
        """ Assembly containing only assigned regions (i.e. regions with contiguous KnobsIntoHoles. """
        g = self.filter_graph(self.graph, cutoff=cutoff, min_kihs=min_kihs)
        ccs = sorted(networkx.connected_component_subgraphs(g, copy=True),
                                 key=lambda x: len(x.nodes()), reverse=True)
        cc = ccs[cc_number]
        helices = [x for x in g.nodes() if x.number in cc.nodes()]
        assigned_regions = self.get_assigned_regions(helices=helices, include_alt_states=False, complementary_only=True)
        coiledcoil_monomers = [h.get_slice_from_res_id(*assigned_regions[h.number]) for h in helices]
        return Assembly(coiledcoil_monomers)

    @property
    def daisy_chain_graph(self):
        """ Directed graph with edges from knob residue to each hole residue for each KnobIntoHole in self. """
        g = networkx.DiGraph()
        for x in self.get_monomers():
            for h in x.hole:
                g.add_edge(x.knob, h)
        return g

    def daisy_chains(self, kih, max_path_length=None):
        """ Generator for daisy chains (complementary kihs) associated with a knob.

        Notes
        -----
        Daisy chain graph is the directed graph with edges from knob residue to each hole residue for each KnobIntoHole
         in self.
        Given a KnobIntoHole, the daisy chains are non-trivial paths in this graph (walks along the directed edges)
        that begin and end at the knob.
        These paths must be of length  <= max_path_length

        Parameters
        ----------
        kih : KnobIntoHole interaction.
        max_path_length : int or None
            Maximum length of a daisy chain.
            Defaults to number of chains in self.ampal_parent.
            This is the maximum sensible value. Larger values than this will cause slow running of this function.
        """
        if max_path_length is None:
            max_path_length = len(self.ampal_parent)
        g = self.daisy_chain_graph
        paths = networkx.all_simple_paths(g, source=kih.knob, target=kih.knob, cutoff=max_path_length)
        return paths

    def get_assigned_regions(self, helices=None, include_alt_states=False, complementary_only=False, cutoff=None):
        """

        Parameters
        ----------
        helices : None, or subset of chains in self.ampal_parent.
            If None, helices is set to self.ampal_parent.
        include_alt_states : bool
            If True, include alternate states for residues.
            If False, include only residues with a single defined state.
        complementary_only : bool
            Consider only KnobIntoHoles with complementary knobs.
        cutoff : None or float
            Cutoff at which to determine complementarity.

        Returns
        -------
        assigned_regions : dict
            Keys : int
                The helix number
            Values : tuple(int)
                The Residue numbers of the start and end of the assigned region.
        """
        if helices is None:
            helices = self.ampal_parent
            kihs = self.get_monomers()
        else:
            kihs = [x for x in self.get_monomers() if (x.knob_helix in helices) and (x.hole_helix in helices)]
        if complementary_only:
            kihs = [x for x in kihs if x.is_complementary(cutoff=cutoff)]
        assigned_regions = {}
        for h in helices:
            kih_residues = set(itertools.chain.from_iterable([x.residues for x in kihs]))
            h_cc_residues = kih_residues.intersection(set(h.get_monomers()))
            if not h_cc_residues:
                continue
            if not include_alt_states:
                h_cc_residues = list(filter(lambda x: len(x.states) == 1, h_cc_residues))
            h_cc_res_numbers = sorted([int(x.id) for x in h_cc_residues])
            assigned_regions[h.number] = (h_cc_res_numbers[0], h_cc_res_numbers[-1])
        return assigned_regions


class KnobIntoHole(PseudoMonomer):

    def __init__(self, pseudo_atoms=None, mol_code='UNK', monomer_id=' ', insertion_code=' ', ampal_parent=None):
        super(KnobIntoHole, self).__init__(pseudo_atoms=pseudo_atoms, monomer_id=monomer_id, ampal_parent=ampal_parent)
        self.mol_code = mol_code
        self.insertion_code = insertion_code
        self.is_hetero = True

    def __repr__(self):
        return '<KnobIntoHole from {0} to [{1}]>'.format(self.knob.name, [x.name for x in self.hole])

    @property
    def knob_helix(self):
        return self.knob.ampal_parent.tags['helix']

    @property
    def hole_helix(self):
        return self.hole[0].ampal_parent.tags['helix']

    @property
    def knob_chain(self):
        return self.knob_helix.ampal_parent.id

    @property
    def hole_chain(self):
        return self.hole_helix.ampal_parent.id

    @property
    def knob_residue(self):
        return self.knob.ampal_parent

    @property
    def hole_residues(self):
        return [x.ampal_parent for x in self.hole]

    @property
    def residues(self):
        return [self.knob_residue] + self.hole_residues

    @property
    def knob(self):
        return self.atoms['k']

    @property
    def hole(self):
        return [v for k, v in self.atoms.items() if 'h' in k]

    @property
    def knob_end(self):
        """ Coordinates of the end of the knob residue (atom in side-chain furthest from CB atom.
        Returns CA coordinates for GLY.
        """
        side_chain_atoms = self.knob_residue.side_chain
        if not side_chain_atoms:
            return self.knob_residue['CA']
        distances = [distance(self.knob_residue['CB'], x) for x in side_chain_atoms]
        max_d = max(distances)
        knob_end_atoms = [atom for atom, d in zip(side_chain_atoms, distances) if d == max_d]
        if len(knob_end_atoms) == 1:
            return knob_end_atoms[0]._vector
        else:
            return numpy.mean([x._vector for x in knob_end_atoms], axis=0)

    @property
    def max_kh_distance(self):
        return max([distance(self.knob, h) for h in self.hole])

    @property
    def packing_angle(self):
        """ Angle between CA-CB of knob and CA(h1)-CA(h2). Returns None if knob is GLY. """
        try:
            knob_vector = self.knob_residue['CB'] - self.knob_residue['CA']
        # exception for GLY residues (with no CB atom).
        except KeyError:
            return None
        hole_vector = self.hole_residues[2]['CA'] - self.hole_residues[1]['CA']
        return angle_between_vectors(knob_vector, hole_vector)

    @property
    def max_knob_end_distance(self):
        """ Maximum distance between knob_end and each of the hole side-chain centres. """
        return max([distance(self.knob_end, h) for h in self.hole])

    @property
    def is_parallel(self):
        v1 = polypeptide_vector(self.knob_helix)
        v2 = polypeptide_vector(self.hole_helix)
        return is_acute(v1, v2)

    def is_complementary(self, cutoff=None):
        knob_group = self.ampal_parent
        if not isinstance(knob_group, KnobGroup):
            raise TypeError("Method only possible when ampal_parent is a KnobGroup. (ampal_parent = {0})"
                            .format(self.ampal_parent))
        if cutoff is not None:
            if cutoff < self.max_kh_distance:
                complementary = False
                return complementary
            else:
                knob_group = knob_group.knob_subgroup(cutoff=cutoff)
        complementary = True if next(knob_group.daisy_chains(kih=self), None) is not None else False
        return complementary

    def knob_type(self, cutoff=None, insertion_cutoff=7.0):
        kt = 1 if self.max_knob_end_distance > insertion_cutoff else 2
        if self.is_complementary(cutoff=cutoff):
            kt += 2
        return kt


def make_pymol(pdb_file, cutoff=7.0, min_kihs=2, outfile=None):
    """ Pymol script for viewing classic coiled-coil Socket output.

    Notes
    -----
    For examples of these views, browse the CC+ database here: http://coiledcoils.chm.bris.ac.uk/ccplus/search/.

    Parameters
    ----------
    pdb_file: str
        Path to a pdb_file.
    cutoff: float
        Socket cutoff in Angstroms.
    min_kihs: int
        Mininmum number of KnobIntoHole interactions between pairs of helices needed to define a coiled coil.
    outfile: None or str
        Path to a output file to save the pml script.

    Returns
    -------
    script_string: str
        Pymol commands for classic coiled-coil view.

    """
    a = convert_pdb_to_ampal(pdb=pdb_file, path=True)
    kg = KnobGroup.from_helices(a, cutoff=cutoff)
    g = kg.filter_graph(kg.graph, cutoff=cutoff, min_kihs=min_kihs)
    ccs = sorted_connected_components(g)
    # Opens pymol script, initial set up of screen
    script_lines = ['load {0}'.format(pdb_file)]
    script_lines.append("hide all")
    script_lines.append("bg_color white")
    script_lines.append("set antialias, 1")
    script_lines.append("set cartoon_dumbbell_length, 0.35")
    script_lines.append("set_color lightgrey, [240,240,240]")
    script_lines.append("set depth_cue, 0")
    script_lines.append("color lightgrey, all")
    script_lines.append("cartoon dumbbell")
    script_lines.append("show cartoon")
    for cc_number, cc in enumerate(ccs):
        helices = [x for x in g.nodes() if x.number in cc.nodes()]
        #helices = cc.nodes()
        cc_region = kg.get_coiledcoil_region(cc_number=cc_number, cutoff=cutoff, min_kihs=min_kihs)
        tag_residues_with_heptad_register(cc_region)
        assigned_regions = kg.get_assigned_regions(include_alt_states=False, complementary_only=False, helices=helices)
        helix_starts = [int(h[0].id) for h in helices]
        helix_ends = [int(h[-1].id) for h in helices]
        chains = [h.ampal_parent.id for h in helices]
        assigned_starts = [assigned_regions[h.number][0] for h in helices]
        assigned_ends = [assigned_regions[h.number][1] for h in helices]
        assigned_selections = ['{0}/{1}-{2}/'.format(chain, assigned_start, assigned_end)
                               for chain, assigned_start, assigned_end in zip(chains, assigned_starts, assigned_ends)]
        script_lines.append("select cc{0}, {1}".format(cc_number, ' '.join(assigned_selections)))
        script_lines.append("cartoon automatic, cc{0}".format(cc_number))
        for h_number, h in enumerate(helices):
            chain = chains[h_number]
            helix_start = helix_starts[h_number]
            helix_end = helix_ends[h_number]
            assigned_start = assigned_starts[h_number]
            assigned_end = assigned_ends[h_number]
            selection = '{0}/{1}-{2}/'.format(chain, helix_start, helix_end)
            script_lines.append("select cc{0}eh{1}, {2}".format(cc_number, h_number, selection))
            selection = '{0}/{1}-{2}/'.format(chain, assigned_start, assigned_end)
            script_lines.append("select cc{0}ah{1}, {2}".format(cc_number, h_number, selection))
            kihs = [x for x in kg if x.knob_helix == h]
            for x in kihs:
                knob_selection_name = 'cc{0}ah{1}k{2}'.format(cc_number, h_number, x.knob_residue.id)
                hole_selection_name = knob_selection_name + 'hole'
                knob_selection = '{0}/{1}/'.format(chain, x.knob_residue.id)
                script_lines.append('select {0}, {1}'.format(knob_selection_name, knob_selection))
                hole_selection = ' '.join(['{0}/{1}/'.format(x.hole_chain, y.id) for y in x.hole_residues])
                script_lines.append('select {0}, {1}'.format(hole_selection_name, hole_selection))
                script_lines.append('show sticks, {0}'.format(knob_selection_name))
                script_lines.append('show sticks, {0}'.format(hole_selection_name))
            for r in h.get_monomers():
                if 'register' in r.tags:
                    color = _heptad_colours[r.tags['register']]
                    script_lines.append('color {0}, {1}/{2}/'.format(color, chain, r.id))
    script_lines.append('deselect')
    script_lines.append('orient')
    script_lines.append('rotate z, 90')
    script_lines.append('zoom complete=1')
    script_string = '\n'.join(script_lines)
    if outfile is not None:
        if isinstance(outfile, str) and outfile[-3:] == 'pml':
            with open(outfile, 'w') as foo:
                foo.write(script_string)
    return script_string


def start_and_end_of_reference_axis(chains):
    """ Get start and end coordinates that approximate the reference axis for a collection of chains
     (not necessarily all the same length).

    Parameters
    ----------
    chains : [Polypeptide]

    Returns
    -------
    start, end : numpy.array
        3D start and end coordinates for defining the reference axis.
    """
    coords = [numpy.array(chains[0].primitive.coordinates)]
    orient_vector = polypeptide_vector(chains[0])
    # Append the coordinates for the remaining chains, reversing the direction in antiparallel arrangements.
    for i, c in enumerate(chains[1:]):
        if is_acute(polypeptide_vector(c), orient_vector):
            coords.append(numpy.array(c.primitive.coordinates))
        else:
            coords.append(numpy.flipud(numpy.array(c.primitive.coordinates)))
    start = numpy.mean([x[0] for x in coords], axis=0)
    end = numpy.mean([x[-1] for x in coords], axis=0)
    return start, end


def gen_reference_primitive(polypeptide, start, end):
    """ Generates a reference Primitive for a Polypeptide given start and end coordinates.

    Notes
    -----
    Uses the rise_per_residue of the Polypeptide primitive to define the separation of points on the line joining
    start and end.

    Parameters
    ----------
    polypeptide : Polypeptide
    start : numpy.array
        3D coordinates of reference axis start
    end : numpy.array
        3D coordinates of reference axis end

    Returns
    -------
    reference_primitive : Primitive
    """
    prim = polypeptide.primitive
    q = find_foot(a=start, b=end, p=prim.coordinates[0])
    ax = Axis(start=q, end=end)
    # flip axis if antiparallel to polypeptide_vector
    if not is_acute(polypeptide_vector(polypeptide), ax.unit_tangent):
        ax = Axis(start=end, end=q)
    arc_length = 0
    points = [ax.start]
    for rise in prim.rise_per_residue()[:-1]:
        arc_length += rise
        t = ax.t_from_arc_length(arc_length=arc_length)
        point = ax.point(t)
        points.append(point)
    reference_primitive = Primitive.from_coordinates(points)
    return reference_primitive


def tag_residues_with_heptad_register(helices):
    """ tags Residues in input helices with heptad register. (Helices not required to be the same length).

    Parameters
    ----------
    helices : [Polypeptide]

    Returns
    -------
    None
    """
    base_reg = 'abcdefg'
    start, end = start_and_end_of_reference_axis(helices)
    for h in helices:
        ref_axis = gen_reference_primitive(h, start=start, end=end)
        crangles = crick_angles(h, reference_axis=ref_axis, tag=False)[:-1]
        reg_fit = fit_heptad_register(crangles)
        exp_base = base_reg * (len(h) // 7 + 2)
        hep_pos = reg_fit[0][0]
        register_string = exp_base[hep_pos:hep_pos + len(h)]
        for i, register in enumerate(register_string):
            h[i].tags['register'] = register
    return


__author__ = 'Jack W. Heal'
