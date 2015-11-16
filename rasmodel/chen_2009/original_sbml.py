import sys
import collections
import re
import inspect
import os
import pkg_resources
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

def load_model():
    """Load model from SBML and return it, also store it in global `model`."""

    global model

    if 'model' in globals():
        return model

    parent_package = __name__.rpartition('.')[0]
    f = pkg_resources.resource_stream(parent_package,
                                      'data/original_sbml_norules.xml')

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


def pysb_model():
    """Return a PySB direct conversion of the SBML model."""

    import pysb
    sbml_model = load_model()

    model = pysb.Model('chen_2009_original_sbml', _export=False)
    for s in sbml_model.species:
        monomer = pysb.Monomer(s.name, _export=False)
        model.add_component(monomer)
        if s.initial_amount != 0:
            parameter = pysb.Parameter(s.name + '_0', s.initial_amount,
                                       _export=False)
            model.add_component(parameter)
            model.initial(monomer, parameter)
    for p in set([p for r in sbml_model.reactions for p in r.kf, r.kr]):
        parameter = pysb.Parameter(p.name, p.value, _export=False)
        model.add_component(parameter)
    for r in sbml_model.reactions:
        assert 1 <= len(r.reactants) <= 2
        assert len(r.products) == 1
        lterms = [model.monomers[s.name]() for s in r.reactants]
        lhs = sum(lterms[1:], lterms[0]())
        rhs = model.monomers[r.products[0].name]()
        kf = model.parameters[r.kf.name]
        kr = model.parameters[r.kr.name]
        if len(lterms) == 2 and r.reactants[0] == r.reactants[1]:
            expr_name = kf.name + '_symmetric'
            try:
                kf = model.expressions[expr_name]
            except KeyError:
                kf = pysb.Expression(expr_name, kf * 2, _export=False)
                model.add_component(kf)
        if kr.get_value() == 0:
            rule = pysb.Rule(r.name, lhs >> rhs, kf, _export=False)
        else:
            rule = pysb.Rule(r.name, lhs <> rhs, kf, kr, _export=False)
        model.add_component(rule)
    sbml_perb11_names = ['c483', 'c136', 'c23', 'c7', 'c25', 'c88', 'c27', 'c89', 'c29', 'c90', 'c34', 'c91', 'c35', 'c92', 'c36', 'c93', 'c37', 'c94', 'c68', 'c67', 'c66', 'c65', 'c21', 'c20', 'c18', 'c19', 'c5', 'c8', 'c15', 'c17', 'c32', 'c63', 'c33', 'c64', 'c95', 'c97', 'c99', 'c419', 'c100', 'c420', 'c486', 'c104', 'c448', 'c415', 'c489', 'c431', 'c264', ]
    sbml_perb12_names = ['c427', 'c130', 'c189', 'c195', 'c198', 'c204', 'c207', 'c213', 'c216', 'c222', 'c225', 'c231', 'c243', 'c249', 'c252', 'c258', 'c234', 'c240', 'c237', 'c255', 'c246', 'c228', 'c219', 'c210', 'c201', 'c192', 'c148', 'c162', 'c165', 'c151', 'c180', 'c183', 'c171', 'c174', 'c445', 'c261', 'c449', 'c416', 'c464', 'c433', 'c265', ]
    sbml_perb13_names = ['c428', 'c131', 'c190', 'c196', 'c199', 'c205', 'c208', 'c214', 'c217', 'c223', 'c226', 'c232', 'c244', 'c250', 'c253', 'c259', 'c235', 'c241', 'c238', 'c256', 'c247', 'c229', 'c220', 'c211', 'c202', 'c193', 'c149', 'c163', 'c166', 'c152', 'c181', 'c184', 'c172', 'c175', 'c446', 'c262', 'c450', 'c281', 'c465', 'c435', 'c409', 'c266', 'c411', ]
    sbml_perb14_names = ['c429', 'c132', 'c191', 'c197', 'c200', 'c206', 'c209', 'c215', 'c218', 'c224', 'c227', 'c233', 'c245', 'c251', 'c254', 'c260', 'c236', 'c242', 'c239', 'c257', 'c248', 'c230', 'c221', 'c212', 'c203', 'c194', 'c150', 'c164', 'c167', 'c153', 'c182', 'c185', 'c173', 'c176', 'c447', 'c263', 'c451', 'c282', 'c466', 'c438', 'c410', 'c267', 'c412', ]
    sbml_perk_names = ['c59', 'c61', 'c95', 'c97', 'c101', 'c431', 'c433', 'c435', 'c438', 'c474', 'c477', 'c480']
    sbml_pakt_names = ['c497', 'c498', 'c472']
    pErbB11_cp = model.monomers[sbml_perb11_names[0]]()
    for n in sbml_perb11_names[1:]:
        pErbB11_cp += model.monomers[n]()
    pErbB1n_cp = model.monomers[sbml_perb12_names[0]]()
    for n in sbml_perb12_names[1:] + sbml_perb13_names + sbml_perb14_names:
        pErbB1n_cp += model.monomers[n]()
    pERK_cp = model.monomers[sbml_perk_names[0]]()
    for n in sbml_perk_names[1:]:
        pERK_cp += model.monomers[n]()
    pAKT_cp = model.monomers[sbml_pakt_names[0]]()
    for n in sbml_pakt_names[1:]:
        pAKT_cp += model.monomers[n]()
    pErbB11 = pysb.Observable('pErbB11', pErbB11_cp, _export=False)
    pErbB1n = pysb.Observable('pErbB1n', pErbB1n_cp, _export=False)
    pErbB1 = pysb.Expression('pErbB1', pErbB11 * 2 + pErbB1n, _export=False)
    pERK = pysb.Observable('pERK', pERK_cp, _export=False)
    pAKT = pysb.Observable('pAKT', pAKT_cp, _export=False)
    for c in pErbB11, pErbB1n, pErbB1, pERK, pAKT:
        model.add_component(c)
    return model
