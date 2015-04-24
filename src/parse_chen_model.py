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

cluster_seq = itertools.count()
def add_cluster(*species_list):
    nodes = [s.id for s in species_list]
    name = 'cluster_{}'.format(next(cluster_seq))
    graph.add_subgraph(nodes, name, color='none', bgcolor='gray97')

graph = pygraphviz.AGraph(directed=True, rankdir='LR')
graph.node_attr.update(fontname='Helvetica')
for s in species:
    graph.add_node(s.id, _type='species', label=s.label, shape='none', bgcolor='white', margin=0)
for r in reactions:
    graph.add_node(r.id, _type='reaction', label=r.name, shape='none', fontcolor='#13ac4a')
    for reactant in r.reactants:
        graph.add_edge(reactant.id, r.id)
    for product in r.products:
        graph.add_edge(r.id, product.id)

# single receptors and 2:3 / 2:4 dimers
add_cluster(c141, c140, c143, c531, c288, c117)
# ligands
add_cluster(c1, c514)
# ligand-bound single receptors
add_cluster(c3, c142, c144)
# ligand-bound receptor dimers
add_cluster(c4, c145, c146, c147, c355, c345, c516, c517)
# ATP-and-ligand-bound receptor dimers
add_cluster(c116, c122, c127, c128, c168, c139, c137, c138)

# delete some "nuisance" nodes from the graph
#for nuisance in ('ATP', 'R_degraded', 'Inh'):
#   graph.remove_node(species_by_label(nuisance).id)

clustered_species = [n for g in graph.subgraphs() for n in g.nodes()]
cs_edges = graph.edges(clustered_species)
nodes_keep = set(reduce(operator.add, map(list, cs_edges), []))
nodes_drop = set(graph.nodes()) - nodes_keep
nd_edges = graph.edges(nodes_drop)
nodes_drop.update(set(reduce(operator.add, map(list, nd_edges), [])))
for n in list(nodes_drop):
    graph.remove_node(n)

graph.write('chen_2009.dot')
