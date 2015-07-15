import sys
import collections
import re
import inspect
import os
import lxml.etree


def simple_repr(x):
    fields = list(x._fields)
    fields.remove('id')
    fields_repr = ', '.join('{0}={1!r}'.format(f, getattr(x, f))
                            for f in fields)
    return '{0}({1})'.format(type(x).__name__, fields_repr)

class Species(collections.namedtuple(
        'SpeciesBase',
        'id name label compartment initial_amount constant')):

    __repr__ = simple_repr

    def __str__(self):
        return self.name

class Parameter(collections.namedtuple(
        'ParameterBase',
        'id name value constant notes')):

    __repr__ = simple_repr

    def __str__(self):
        return self.name

class Reaction(collections.namedtuple(
        'ReactionBase',
        'id name label reactants products kf kr')):

    def __str__(self):
        return self.name

    def __repr__(self):
        simple_attrs = ('name', 'label')
        simple_attr_repr = ', '.join('{0}={1!r}'.format(a, getattr(self, a))
                                     for a in simple_attrs)
        species_repr = ', '.join(
            ['{0}=({1})'.format(et, ', '.join(s.name for s in getattr(self, et)))
             for et in 'reactants', 'products'])
        param_repr = ', '.join(['{0}={1}'.format(pt, getattr(self, pt))
                                for pt in 'kf', 'kr'])
        return '{0}({1}, {2}, {3})'.format(type(self).__name__,
                                           simple_attr_repr, species_repr,
                                           param_repr)

class Model(collections.namedtuple(
        'ModelBase',
        'species parameters reactions')):

    def species_by_label(self, pattern, multi=False):
        if isinstance(pattern, basestring):
            pattern = r'^{}$'.format(re.escape(pattern))
        gen = (s for s in self.species if re.search(pattern, s.label))
        return list(gen) if multi else next(gen)

    def export_globals(self, dest=None):
        if dest is None:
            dest = inspect.currentframe().f_back.f_globals
        # First collect elements in a dict rather than updating the dest
        # namespace directly. This lets us raise an exception if there's a
        # duplicate without leaving the dest namespace in a partially updated
        # state.
        elements = {}
        for container in self:
            for elt in container:
                if elt.name in elements:
                    raise RuntimeError('duplicate component name: {}'
                                       .format(elt.name))
                elements[elt.name] = elt
        # Now that we're sure there are no duplicates, update the namespace in
        # one shot.
        dest.update(elements)

    def get_duplicate_names(self):
        names = set()
        duplicates = set()
        for container in self:
            for elt in container:
                if elt.name in names:
                    duplicates.add(elt.name)
                else:
                    names.add(elt.name)
        return duplicates


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

def reaction_species(element, entity_type):
    entity_type = entity_type.title()
    query = 's:listOf{}s/s:speciesReference/@species'.format(entity_type)
    return xpath(element, query, single=False)


ns = 'http://www.sbml.org/sbml/level2'
qnames = dict((tag, lxml.etree.QName(ns, tag).text)
              for tag in ('species', 'reaction', 'parameter'))

label_edits = {
    'c95': '2(EGF:ErbB1)#P:GAP:Grb2:Sos:ERK#P#P',
    'c96': 'endo|2(EGF:ErbB1)#P:GAP:Grb2:Sos:ERK#P#P',
    'c98': 'endo|2(EGF:ErbB1)#P:GAP:(Shc#P):Grb2:Sos:ERK#P#P',
    'c100': 'endo|2(EGF:ErbB1)#P:GAP:Grb2:Sos#P',
    'c101': 'Sos:ERK#P#P',
    'c105': 'ATP',
    'c123': 'endo|EGF:ErbB1:ErbB2:ATP',
    'c124': 'endo|EGF:ErbB1:ErbB3:ATP',
    'c125': 'endo|EGF:ErbB1:ErbB4:ATP',
    'c126': 'endo|2(EGF:ErbB1:ATP)-FullActive',
    'c157': 'endo|HRG:ErbB3',
    'c158': 'endo|HRG:ErbB4',
    'c169': 'endo|(HRG:ErbB3):ErbB2:ATP',
    'c170': 'endo|(HRG:ErbB4):ErbB2:ATP',
    'c264': '2(EGF:ErbB1)#P:GAP:Grb2:Gab1#P:PI3K:Ras:GDP',
    'c284': 'ErbB2#P:ErbB2',
    'c288': 'ErbB2:ErbB3',
    'c339': 'endo|ErbB2:ErbB3',
    'c340': 'endo|ErbB2:ErbB4',
    'c407': '(ErbB4:ErbB2)#P:GAP:Grb2:Gab1#P#P:Pase9t',
    'c409': '(ErbB1:ErbB3)#P:GAP:Grb2:Gab1#P#P',
    'c410': '(ErbB1:ErbB4)#P:GAP:Grb2:Gab1#P#P',
    'c411': '(ErbB1:ErbB3)#P:GAP:Grb2:Gab1#P#P:Pase9t',
    'c412': '(ErbB1:ErbB4)#P:GAP:Grb2:Gab1#P#P:Pase9t',
    'c417': 'endo|(ErbB3:ErbB2)#P:RTK_Pase',
    'c418': 'endo|(ErbB4:ErbB2)#P:RTK_Pase',
    'c419': '2(EGF:ErbB1)#P:GAP:(Shc#P):Grb2:Sos#P',
    'c420': 'endo|2(EGF:ErbB1)#P:GAP:(Shc#P):Grb2:Sos#P',
    'c421': 'endo|(HRG:ErbB3):ErbB2',
    'c422': 'endo|(HRG:ErbB4):ErbB2',
    'c424': '(ErbB1:ErbB2)#P:GAP:Grb2:Gab1#P#P:Pase9t',
    'c425': 'endo|ErbB2:ErbB2',
    'c430': '(ErbB1:ErbB2)#P:GAP:Grb2:Gab1#P#P',
    'c431': '2(EGF:ErbB1)#P:GAP:Grb2:Gab1#P:ERK#P#P',
    'c438': '(ErbB1:ErbB4)#P:GAP:Grb2:Gab1#P:ERK#P#P',
    'c455': '(ErbB4:ErbB2)#P:GAP:Grb2:Gab1#P:PI3K:PIP2',
    'c456': '(ErbB3:ErbB2)#P:GAP:Grb2:Gab1#P#P:Pase9t',
    'c472': 'Raf#P:AKT:P:P',
    'c486': '2(EGF:ErbB1)#P:GAP:Grb2:Gab1#P',
    'c487': '(ErbB4:ErbB2)#P:GAP:Grb2:Gab1#P#P',
    'c488': '2(EGF:ErbB1)#P:GAP:Grb2:Gab1#P#P',
    'c489': '2(EGF:ErbB1)#P:GAP:Grb2:Gab1#P:Shp2',
    'c490': '2(ErbB2)#P:GAP:Grb2:Gab1#P#P',
    'c491': '(ErbB3:ErbB2)#P:GAP:Grb2:Gab1#P#P',
    'c522': '2(EGF:ErbB1)#P:GAP:Grb2:Gab1#P#P:Pase9t',
    'c523': '2(ErbB2)#P:GAP:Grb2:Gab1#P#P:Pase9t',
    }

def get_model():

    model_filename = os.path.join(os.path.dirname(__file__),
                                  'chen_2009_original_sbml_norules.xml')
    f = open(model_filename)

    model = Model([], [], [])
    species_id_map = {}
    parameters_name_map = {}

    for event, element in lxml.etree.iterparse(f, tag=qnames['species']):
        species_id = element.get('id')
        name = element.get('name')
        label = xpath(element, 's:notes/text()').strip()
        compartment = xpath(element, 's:annotation/text()')
        if compartment is not None:
            compartment = compartment.strip().lower()
            if 'endo' in compartment:
                label = 'endo|' + label
        initial_amount = float(element.get('initialAmount'))
        const = element.get('constant')
        if const == 'true':
            const = True
        elif const is None:
            const = False
        else:
            raise RuntimeError('bad species.constant value: {}'.format(const))
        # Override label for some badly named species.
        if name in label_edits:
            label = label_edits[name]
        s = Species(species_id, name, label, compartment, initial_amount, const)
        model.species.append(s)
        species_id_map[s.id] = s

    f.seek(0)
    for event, element in lxml.etree.iterparse(f, tag=qnames['parameter']):
        parameter_id = element.get('id')
        name = element.get('name')
        value = float(element.get('value'))
        const = element.get('constant')
        if const is None:
            const = True
        elif const == 'false':
            const = False
        else:
            raise RuntimeError('bad parameter.constant value: {}'.format(const))
        notes = xpath(element, 's:notes/text()')
        if notes is not None:
            notes = notes.strip()
        p = Parameter(parameter_id, name, value, const, notes)
        model.parameters.append(p)
        parameters_name_map[p.name] = p

    f.seek(0)
    for event, element in lxml.etree.iterparse(f, tag=qnames['reaction']):
        rxn_id = element.get('id')
        label = re.sub(r' +', ' ', element.get('name'))
        name, label = label.split(' ', 1)
        label, kf, kr = label.rsplit(' ', 2)
        kf = parameters_name_map[kf]
        kr = parameters_name_map[kr]
        reactants = tuple(map(species_id_map.get,
                              reaction_species(element, 'reactant')))
        products = tuple(map(species_id_map.get,
                             reaction_species(element, 'product')))
        r = Reaction(rxn_id, name, label, reactants, products, kf, kr)
        model.reactions.append(r)

    duplicates = model.get_duplicate_names()
    if duplicates:
        raise RuntimeError('duplicate names: ' + ', '.join(duplicates))

    return model
