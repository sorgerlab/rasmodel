from __future__ import division
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

def reverse_reaction(reaction, keep_appearance=False):
    for s in reaction.reactants:
        if s.id in graph:
            reverse_edge(s.id, reaction.id, keep_appearance)
    for s in reaction.products:
        if s.id in graph:
            reverse_edge(reaction.id, s.id, keep_appearance)

def reverse_edge(u, v, keep_appearance):
    e = graph.get_edge(u, v)
    attrs = dict(e.attr)
    if keep_appearance:
        attrs['dir'] = 'back'
    graph.add_edge(v, u, **attrs)
    graph.remove_edge(u, v)

cluster_seq = itertools.count()
def add_box(*species_list):
    nodes = [s.id for s in species_list]
    name = 'cluster_{}'.format(next(cluster_seq))
    cluster = graph.add_subgraph(nodes, name, color='none', bgcolor='gray94')
    return cluster

# All the dark ("4") colors from the graphviz color set.
color_cycle = itertools.cycle((
        'antiquewhite4', 'aquamarine4', 'azure4', 'bisque4', 'blue4', 'brown4',
        'burlywood4', 'cadetblue4', 'chartreuse4', 'chocolate4', 'coral4',
        'cornsilk4', 'cyan4', 'darkgoldenrod4', 'darkolivegreen4',
        'darkorange4', 'darkorchid4', 'darkseagreen4', 'deeppink4',
        'deepskyblue4', 'dodgerblue4', 'firebrick4', 'gold4', 'goldenrod4',
        'green4', 'honeydew4', 'hotpink4', 'indianred4', 'ivory4', 'khaki4',
        'lavenderblush4', 'lemonchiffon4', 'magenta4', 'maroon4',
        'mediumorchid4', 'mediumpurple4', 'mistyrose4', 'navajowhite4',
        'olivedrab4', 'orange4', 'orangered4', 'orchid4', 'peachpuff4', 'pink4',
        'plum4', 'purple4', 'red4', 'rosybrown4', 'royalblue4', 'salmon4',
        'seagreen4', 'seashell4', 'sienna4', 'skyblue4', 'slateblue4', 'snow4',
        'springgreen4', 'steelblue4', 'tan4', 'thistle4', 'tomato4',
        'turquoise4', 'violetred4', 'wheat4', 'yellow4'))

#####
# Parse SBML file.

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
        if 'endo' in compartment:
            label = 'endo|' + label
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

#####
# Construct graphviz representation of reaction network.

graph = pygraphviz.AGraph(directed=True, rankdir='LR', compound=True)
graph.node_attr.update(fontname='Helvetica')
graph.edge_attr.update(arrowsize=0.7)
for s in species:
    label = '<{0.label} <sup>{0.name}</sup>>'.format(s)
    graph.add_node(s.id, _type='species', label=label, shape='none',
                   bgcolor='white', margin=0.01, height=0)
for r in reactions:
    graph.add_node(r.id, _type='reaction', label=r.name, shape='none',
                   fontcolor='#13ac4a', width=0, height=0, margin=0.05)
    color = next(color_cycle)
    for reactant in r.reactants:
        graph.add_edge(reactant.id, r.id, color=color)
    for product in r.products:
        graph.add_edge(r.id, product.id, color=color)

# delete some "nuisance" nodes from the graph
for label in 'ATP',:
    graph.remove_node(species_by_label(label).id)

# Fix reactions that were specified "backwards" in the original model. This
# swaps the direction of all edges on these reaction nodes.
for r in (
    v804, v807, v812, v813, v822, v823, v824, v825, # dimer#P catalysis step
    v808, v811, v809, v826, v810, v827,       # endo|dimer#P catalysis step
    v657, v658, v659, v660, v661, v662, v663, # endo|dimer#P diss RTK_Pase
    v109, v111, v123, v139, v140, v141, v161, # R.. diss cPP
    v107, v108, v122, v128, v129, v130, v162, # R..Sos diss cPP
    v117, v118, v127, v151, v152, v153, v158, # R..Shc diss cPP
    v110, v116, v126, v148, v149, v150, v157, # R..Shc..Sos diss cPP
    v105, v106, v121, v136, v137, v138, v160, # R..RasGDP diss cPP
    v103, v104, v120, v133, v134, v135, v159, # R..RasGTP diss cPP
    v114, v115, v125, v145, v146, v147, v156, # R..Shc..RasGDP diss cPP
    v112, v113, v124, v142, v143, v144, v155, # R..Shc..RasGTP diss cPP
    v289, v293, v295, v296, v297, v307, v309, # R..Shc..RasGDP diss RasGTP
    v283, v285, v287, v294, v301, v302, v303, # R..RasGDP diss RasGTP
    ):
    reverse_reaction(r)
# Fix "retrograde" reactions -- those whose forward direction is logically
# "backward" in the signaling network. This switches the logical edge direction
# without changing the visual appearance (i.e. which end has the arrowhead),
# leading to improved graph layout.
for r in (v443, v211):
    reverse_reaction(r, keep_appearance=True)

## Plasma membrane receptors and ligands

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

## Adapters

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

## Endosome receptors and ligands

# single receptors
add_box(c156, c154, c155, c6)
# 2:2 / 2:3 / 2:4 dimers
add_box(c425, c339, c340)
# ligands
add_box(c16, c515)
# degraded EGF
add_box(c13)
# ligand-bound single receptors
add_box(c10, c157, c158)
# ligand-bound receptor dimers
add_box(c11, c159, c160, c161, c421, c422, c518, c519)
# ATP-and-ligand-bound receptor dimers
add_box(c126, c123, c124, c125, c169, c170)
# Phosphorylated dimers (dimer#P)
add_box(c8, c162, c163, c164, c337, c338, c290)

## RTK phosphatase

# RTK Phosphatase
add_box(c280)
# dimer#P:Phosphatase
add_box(c415, c416, c281, c283, c282, c417, c418)

## Endosome adapters

# dimer#P:GAP
add_box(c343, c346, c293, c17, c165, c166, c167)
# dimer#P:GAP:Shc
add_box(c349, c350, c296, c174, c175, c176, c63)
# dimer#P:GAP:(Shc#P)
add_box(c353, c356, c299, c183, c184, c185, c64)
# dimer#P:GAP:(Shc#P):Grb2
add_box(c359, c362, c302, c192, c193, c194, c65)
# dimer#P:GAP:(Shc#P):Grb2:Sos
add_box(c365, c368, c305, c201, c202, c203, c66)
# dimer#P:GAP:Grb2
add_box(c383, c386, c314, c228, c229, c230, c18)
# dimer:P:GAP:Grb2:Sos
add_box(c389, c392, c317, c237, c238, c239, c19)

## Coated-pit-protein-bound adapters

# cPP
add_box(c12)
# endo|cPP
add_box(c9)
# dimer#P:GAP:(Shc#P):Grb2
add_box(c358, c361, c301, c195, c196, c197, c91)
# dimer#P:GAP:(Shc#P):Grb2:Sos
add_box(c364, c367, c304, c204, c205, c206, c92)
# dimer#P:GAP:Grb2
add_box(c382, c385, c313, c231, c232, c233, c7)
# dimer:P:GAP:Grb2:Sos
add_box(c388, c391, c316, c240, c241, c242, c88)

## Inhibitor

# Inh
add_box(c285)
# receptor:Inh
add_box(c286, c502, c506, c503)

## RAS-bound complexes

# Ras:GDP
add_box(c26)
# Ras:GTP
add_box(c28)

# dimer#P:GAP:(Shc#P):Grb2:Sos:RasGDP, plasma membrane
add_box(c369, c372, c306, c207, c208, c209, c36)
# dimer#P:GAP:(Shc#P):Grb2:Sos:RasGDP, cpp-bound
add_box(c370, c373, c307, c213, c214, c215, c93)
# dimer#P:GAP:(Shc#P):Grb2:Sos:RasGDP, endosome
add_box(c371, c374, c308, c210, c211, c212, c67)

# dimer#P:GAP:(Shc#P):Grb2:Sos:RasGTP, plasma membrane
add_box(c375, c378, c309, c216, c217, c218, c37)
# dimer#P:GAP:(Shc#P):Grb2:Sos:RasGTP, cpp-bound
add_box(c376, c379, c310, c222, c223, c224, c94)
# dimer#P:GAP:(Shc#P):Grb2:Sos:RasGTP, endosome
add_box(c377, c380, c311, c219, c220, c221, c68)

# dimer#P:GAP:Grb2:Sos:RasGDP, plasma membrane
add_box(c393, c396, c318, c243, c244, c245, c27)
# dimer#P:GAP:Grb2:Sos:RasGDP, cpp-bound
add_box(c394, c397, c319, c249, c250, c251, c89)
# dimer#P:GAP:Grb2:Sos:RasGDP, endosome
add_box(c395, c398, c320, c246, c247, c248, c20)

# dimer#P:GAP:Grb2:Sos:RasGTP, plasma membrane
add_box(c399, c402, c321, c252, c253, c254, c29)
# dimer#P:GAP:Grb2:Sos:RasGTP, cpp-bound
add_box(c400, c403, c322, c258, c259, c260, c90)
# dimer#P:GAP:Grb2:Sos:RasGTP, endosome
add_box(c401, c404, c323, c255, c256, c257, c21)

# Delete stuff we haven't explicitly enumerated through add_box calls above.
box_nodes = [n for g in graph.subgraphs() for n in g.nodes()]
nodes_keep = set(box_nodes).union(neighbor_set(box_nodes))
nodes_drop = set(graph.nodes()) - nodes_keep
dropped_species_nodes = set(n for n in nodes_drop if n.attr['_type'] == 'species')
drop2 = neighbor_set(dropped_species_nodes)
for n in list(nodes_drop.union(drop2)):
    graph.remove_node(n)
num_keep_species = len([n for n in graph.nodes() if n.attr['_type'] == 'species'])
num_keep_reactions = len([n for n in graph.nodes() if n.attr['_type'] == 'reaction'])

# Collapse sets of parallel reaction nodes.
node_to_subgraph = {n: g for g in graph.subgraphs() for n in g.nodes()}
rxn_nodes = [n for n in graph.nodes() if n.attr['_type'] == 'reaction']
rxn_cluster_neighbors = [
    tuple(sorted(node_to_subgraph[sn] for sn in graph.neighbors(rn)))
    for rn in rxn_nodes]
rxn_to_cn = dict(zip(rxn_nodes, rxn_cluster_neighbors))
# Sort reaction nodes based on which species subgraphs (clusters) they are
# adjacent to. Map subgraphs to their ids because graph equality in pygraphviz
# is VERY expensive -- it serializes them to strings and compares those!
rxn_nodes.sort(key=lambda r: map(id, rxn_to_cn[r]))
for subgraphs, node_iter in itertools.groupby(rxn_nodes, rxn_to_cn.get):
    nodes = list(node_iter)
    if (all(sum(n.attr['_type'] == 'species' for n in sg) in (1, len(nodes)) for sg in subgraphs) and
        len(nodes) > 1):
        node_id = 'reaction_' + '_'.join(sg.name for sg in subgraphs)
        label = r'\n'.join(sorted(n.attr['label'] for n in nodes))
        graph.add_node(node_id, label=label, fontcolor='#13ac4a', shape='box',
                       width=0, height=0, margin=0.05, color='#13ac4a20')
        r_nodes = graph.predecessors(nodes[0])
        p_nodes = graph.successors(nodes[0])
        base_edge_attrs = {'color': next(color_cycle)}
        for ntype, species_nodes in ('reactant', r_nodes), ('product', p_nodes):
            for species_node in species_nodes:
                sg = node_to_subgraph[species_node]
                u, v = sg.nodes_iter().next(), node_id
                if ntype == 'reactant':
                    lheadtail = 'ltail'
                else:
                    lheadtail = 'lhead'
                    u, v = v, u
                attrs = dict(base_edge_attrs)
                if len(sg) > 1:
                    attrs.update({lheadtail: sg.name}, arrowsize=1.4,
                                 style='bold')
                graph.add_edge(u, v, **attrs)
        for n in nodes:
            graph.remove_node(n)

# TEST
#graph.add_subgraph([s.id for s in species if s.id in graph and 'endo' in s.compartment],
#                   'endosome', color='red')

# TEMP
# for sg in graph.subgraphs():
#     sg.graph_attr['label'] = sg.name
#     globals()[sg.name] = sg

graph.write('chen_2009.dot')

num_total_species = len(species)
num_total_reactions = len(reactions)
print >>sys.stderr, ("Species: {} / {} ({}%)"
                     .format(num_keep_species, num_total_species,
                             100 * num_keep_species / num_total_species))
print >>sys.stderr, ("Reactions: {} / {} ({}%)"
                     .format(num_keep_reactions, num_total_reactions,
                             100 * num_keep_reactions / num_total_reactions))

#####
# Debugging/cleanup checks

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
