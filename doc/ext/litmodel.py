# Sphinx extension to enable "literate modeling" features for REM.

import sys
import os
import sphinx.addnodes as addnodes
from sphinx.directives.other import TocTree
from sphinx.builders import Builder
from sphinx.util.nodes import inline_all_toctrees
from sphinx.writers.text import TextWriter
from sphinx.writers.html import HTMLWriter
from sphinx.util.console import bold, darkgreen
from sphinx.util.osutil import ensuredir
from sphinx.locale import _
from docutils.parsers.rst import Directive
from docutils.parsers.rst.directives import unchanged_required
from docutils.nodes import note, paragraph, literal, Text, GenericNodeVisitor
from docutils.transforms import Transform
from docutils.writers import Writer

def setup(app):
    app.add_directive('component', ComponentDirective)
    app.add_transform(LineNumberKeeper)
    app.add_builder(TangleBuilder)


class ComponentDirective(Directive):
    """New directive to flag a file as the root of a model component.

    The one required option, 'module', specifies the Python module to which the
    component's source code will be extracted, as a dot-delimited path.

    The directive will be replaced by an informative note in the output document
    which states the module name.
    """

    option_spec = {'module': unchanged_required}

    def run(self):
        module_name = self.options['module']

        # FIXME Check for nested modules? Currently caught in
        # ComponentTranslator below but we should probably do it here.

        env = self.state.document.settings.env
        if not hasattr(env, 'litmodel_components'):
            env.litmodel_components = {}
        env.litmodel_components[env.docname] = module_name

        note_node = note()
        paragraph_node = paragraph('', '',
                                   Text('This section defines the component '),
                                   literal('', module_name),
                                   Text('.'))
        note_node += paragraph_node

        return [note_node]


class LineNumberKeeper(Transform):
    """Preserve node.line values across copy/deepcopy calls.

    node.copy creates the new returned node by calling the appropriate class
    constructor and propagating the 'rawsource' and 'attributes' properties,
    however the 'line' property is not reinitialized. This transform stashes the
    'line' property of every node in a new attribute 'orig_line', and since
    'attributes' is preserved across copies, we can access this attribute
    instead of the line property."""

    default_priority = 400

    def apply(self):
        for node in self.document.traverse():
            if node.line is not None:
                node['orig_line'] = node.line


class TangleBuilder(Builder):
    """Builder that exports components to actual .py files for import."""

    name = 'tangle'
    format = 'tangle'
    out_suffix = '.py'

    def get_outdated_docs(self):
        return 'all documents'

    def get_target_uri(self, docname, typ=None):
        # FIXME Do we need to actually return something useful here or is this OK?
        return ''

    def assemble_doctree(self):
        master = self.config.master_doc
        tree = self.env.get_doctree(master)
        # Final arg is the color function. There is no built-in "no change"
        # color function so we use a trivial lambda expression instead.
        tree = inline_all_toctrees(self, set(), master, tree, lambda x: x)
        tree['docname'] = master
        self.env.resolve_references(tree, master, self)
        return tree

    def component_sof_filter(self, node):
        return (isinstance(node, addnodes.start_of_file) and
                node['docname'] in self.env.litmodel_components)

    def component_sof_stringify(self, node):
        return self.env.litmodel_components[node['docname']]

    def write(self, *ignored):
        # Build a single combined doctree by inlining the toctrees, then iterate
        # over the start_of_file nodes corresponding to the component directives.        

        self.info(bold('preparing documents... '), nonl=True)
        self.prepare_writing()
        self.info('done')

        self.info(bold('assembling single document... '), nonl=True)
        doctree = self.assemble_doctree()
        self.info()

        sof_iterator = doctree.traverse(self.component_sof_filter)
        si = self.status_iterator(sof_iterator, 'writing components... ',
                                  darkgreen, len(self.env.litmodel_components),
                                  self.component_sof_stringify)
        for sof_node in si:
            module_name = self.env.litmodel_components[sof_node['docname']]
            # Writer infrastructure expects a proper document. We shallow-
            # copy our master document, then deep-copy our subtree of
            # interest as its children.
            subdocument = doctree.copy()
            subdocument.children = [sof_node.deepcopy()]
            self.write_module(module_name, subdocument)

    def prepare_writing(self):
        self.docwriter = ComponentWriter(self)

    def write_module(self, name, doctree):
        module_rel_path = name.replace('.', '/')
        # Note that we don't create __init__.py files here, so we can't output
        # anything but "leaf" modules.
        path = os.path.join(self.outdir, module_rel_path) + '.py'
        path = os.path.abspath(path)
        ensuredir(os.path.dirname(path))
        with open(path, 'w') as destination:
            self.docwriter.write(doctree, destination)


class ComponentWriter(Writer):
    # Just enough to bridge to the ComponentTranslator.

    def __init__(self, builder):
        Writer.__init__(self)
        self.builder = builder
        self.translator_class = ComponentTranslator

    def translate(self):
        visitor = self.translator_class(self.document, self.builder)
        self.document.walkabout(visitor)
        self.output = visitor.get_content()


class ComponentTranslator(GenericNodeVisitor):
    """Extracts and merges all literal_blocks."""

    def __init__(self, document, builder):
        GenericNodeVisitor.__init__(self, document)
        self.builder = builder
        self.content = []
        self.docs_with_visited_code = set()
        self.nesting_level = 1

    def depart_start_of_file(self, node):
        if node['docname'] in self.docs_with_visited_code:
            self.nesting_level -= 1
            nesting = '<' * self.nesting_level
            comment = "# %s END %s.rst" % (nesting, node['docname'])
            #self.content.append(comment)

    def visit_literal_block(self, node):
        docname = None
        walk_node = node
        while docname is None and walk_node.parent is not None:
            docname = walk_node.get('docname')
            walk_node = walk_node.parent
        if docname not in self.docs_with_visited_code:
            nesting = '>' * self.nesting_level
            comment = "\n# %s BEGIN %s.rst" % (nesting, docname)
            #self.content.append(comment)
            self.docs_with_visited_code.add(docname)
            self.nesting_level += 1
        self.content.append("\n# SOURCE: %s.rst:%d" % (docname, node['orig_line']))
        self.content.append(node.astext())

    def default_visit(self, node):
        pass

    def default_departure(self, node):
        pass

    def get_content(self):
        return '\n'.join(self.content) + '\n'
