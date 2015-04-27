import sys
import collections
import re
import itertools
import operator
import lxml.etree
import pygraphviz

Species = collections.namedtuple('Species',
                                 'id name label compartment')

class Reaction(collections.namedtuple(
        'ReactionBase',
        'id name label reactants products kf kr')):
    def __repr__(self):
        simple_attrs = ('id', 'name', 'label')
        simple_attr_repr = ', '.join('{0}={1!r}'.format(a, getattr(self, a))
                                     for a in simple_attrs)
        species_repr = ', '.join(
            ['{0}=({1})'.format(et, ', '.join(s.name for s in getattr(self, et)))
             for et in 'reactants', 'products'])
        param_repr = ', '.join(['{0}={1!r}'.format(pt, getattr(self, pt))
                                for pt in 'kf', 'kr'])
        return '{0}({1}, {2}, {3})'.format(type(self).__name__,
                                           simple_attr_repr, species_repr,
                                           param_repr)

def xpath(obj, path, single=True):
    result = map(unicode, obj.xpath(path, namespaces={'s': ns}))
    if single:
        if len(result) == 0:
            result = None
        elif len(result) == 1:
            result = result[0]
        else:
            raise ValueError("XPath expression matches more than one value")
    return result

def reactionSpecies(element, entity_type):
    entity_type = entity_type.title()
    query = 's:listOf{}s/s:speciesReference/@species'.format(entity_type)
    return xpath(element, query, single=False)

def species_by_label(pattern, multi=False):
    if isinstance(pattern, basestring):
        pattern = r'^{}$'.format(re.escape(pattern))
    gen = (s for s in species if re.search(pattern, s.label))
    return list(gen) if multi else next(gen)

def neighbor_set(nodes):
    return set(itertools.chain.from_iterable(graph.neighbors(n) for n in nodes))

def reverse_edge(u, v):
    e = graph.get_edge(u, v)
    graph.add_edge(v, u, **e.attr)
    graph.remove_edge(u, v)

cluster_seq = itertools.count()
def add_box(*species_list):
    nodes = [s.id for s in species_list]
    # failed attempt to force ordering within clusters
    # for n1, n2 in zip(nodes[:-1], nodes[1:]):
    #     graph.add_edge(n1, n2, color='red')
    name = 'cluster_{}'.format(next(cluster_seq))
    return graph.add_subgraph(nodes, name, color='none', bgcolor='gray94')

ns = 'http://www.sbml.org/sbml/level2'
qnames = dict((tag, lxml.etree.QName(ns, tag).text)
              for tag in ('species', 'reaction'))

sbml_file = open(sys.argv[1])

species = []
reactions = []
species_id_map = {}
for event, element in lxml.etree.iterparse(sbml_file, tag=qnames['species']):
    species_id = element.get('id')
    name = element.get('name')
    label = xpath(element, 's:notes/text()').strip()
    # special fixup for weird ATP label
    if label == 'ATP  1.2e9':
        label = 'ATP'
    compartment = xpath(element, 's:annotation/text()')
    if compartment is not None:
        compartment = compartment.strip().lower()
    s = Species(species_id, name, label, compartment)
    species.append(s)
    species_id_map[s.id] = s
    if s.name in globals():
        raise RuntimeError('duplicate component name: {}'.format(s.name))
    globals()[s.name] = s
sbml_file.seek(0)
for event, element in lxml.etree.iterparse(sbml_file, tag=qnames['reaction']):
    rxn_id = element.get('id')
    label = re.sub(r' +', ' ', element.get('name'))
    name, label = label.split(' ', 1)
    label, kf, kr = label.rsplit(' ', 2)
    reactants = tuple(map(species_id_map.get, reactionSpecies(element, 'reactant')))
    products = tuple(map(species_id_map.get, reactionSpecies(element, 'product')))
    r = Reaction(rxn_id, name, label, reactants, products, kf, kr)
    reactions.append(r)
    if r.name in globals():
        raise RuntimeError('duplicate component name: {}'.format(r.name))
    globals()[r.name] = r

# determine truly redundant species -- same name, same compartment
sa = sorted(species, key=lambda s: (s.label, s.compartment))
rs = [x
      for x in ((n, map(lambda s: s.compartment, it))
          for n, it in itertools.groupby(sa, lambda s: s.label))
      if len(x[1]) == 2 and x[1][0] == x[1][1]]

# find reactions where product is not a trivial concatenation of reactants
# (e.g. A + B -> A:B)
mismatch_rxns = [r for r in reactions
                 if r.products[0].label != ':'.join([s.label for s in r.reactants])]

graph = pygraphviz.AGraph(directed=True, rankdir='LR')
graph.node_attr.update(fontname='Helvetica')
for s in species:
    label = '<{0.label} <sup>{0.name}</sup>>'.format(s)
    graph.add_node(s.id, _type='species', label=label, shape='none',
                   bgcolor='white', margin=0.1)
for r in reactions:
    graph.add_node(r.id, _type='reaction', label=r.name, shape='none',
                   fontcolor='#13ac4a', width=0, height=0, margin=0.05)
    for reactant in r.reactants:
        graph.add_edge(reactant.id, r.id, arrowsize=0.7)
    for product in r.products:
        graph.add_edge(r.id, product.id, arrowsize=0.7)

# delete some "nuisance" nodes from the graph
for label in 'ATP', 'R_degraded', 'Inh':
   graph.remove_node(species_by_label(label).id)

# reverse edges for some reactions which were specified "backwards"
for r in v804, v807, v812, v813, v822, v823, v824, v825:
    for s in r.reactants:
        if s.id in graph:
            reverse_edge(s.id, r.id)
    for s in r.products:
        if s.id in graph:
            reverse_edge(r.id, s.id)

# single ErbB1
add_box(c531)
# single receptors
add_box(c2, c141, c140, c143)
# 2:3 / 2:4 dimers
add_box(c288, c117)
# ligands
add_box(c1, c514)
# ligand-bound single receptors
add_box(c3, c142, c144)
# ligand-bound receptor dimers
add_box(c4, c145, c146, c147, c355, c345, c516, c517)
# ATP-and-ligand-bound receptor dimers
add_box(c116, c122, c127, c128, c168, c139, c137, c138)
# Phosphorylated dimers (dimer#P)
add_box(c5, c148, c149, c150, c335, c336, c289)
# Phosphorylated monomers
add_box(c330, c87, c331, c332)
# GAP
add_box(c14)
# dimer#P:GAP
add_box(c341, c344, c291, c15, c151, c152, c153)
# Shc
add_box(c31)
# Grb2
add_box(c22)
# Sos
add_box(c24)
# Shc#P
add_box(c40)
# (Shc#P):Grb2
add_box(c39)
# (Shc#P):Grb2:Sos
add_box(c38)
# Grb2:SOS
add_box(c30)
# dimer#P:GAP:Shc
add_box(c347, c348, c294, c171, c172, c173, c32)
# dimer#P:GAP:(Shc#P)
add_box(c351, c354, c297, c180, c181, c182, c33)
# dimer#P:GAP:(Shc#P):Grb2
add_box(c357, c360, c300, c189, c190, c191, c34)
# dimer#P:GAP:(Shc#P):Grb2:Sos
add_box(c363, c366, c303, c198, c199, c200, c35)
# dimer#P:GAP:Grb2
add_box(c381, c384, c312, c225, c226, c227, c23)
# dimer:P:GAP:Grb2:Sos
add_box(c387, c390, c315, c234, c235, c236, c25)

box_nodes = [n for g in graph.subgraphs() for n in g.nodes()]
nodes_keep = set(box_nodes).union(neighbor_set(box_nodes))
nodes_drop = set(graph.nodes()) - nodes_keep
dropped_species_nodes = set(n for n in nodes_drop if n.attr['_type'] == 'species')
drop2 = neighbor_set(dropped_species_nodes)
for n in list(nodes_drop.union(drop2)):
    graph.remove_node(n)

graph.write('chen_2009.dot')
