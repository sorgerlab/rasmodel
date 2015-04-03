import sys
import docutils.examples

def extract_literals(node):
    if node.tagname == 'literal_block':
        text = node.astext()
    else:
        text = ''
        for child in node.children:
            text += extract_literals(child)
    if text:
        text += '\n\n'
    return text

with open(sys.argv[1]) as f:
    input_string = f.read().decode('utf-8')
settings = {'report_level': 5}
doc, publisher = docutils.examples.internals(input_string,
                                             settings_overrides=settings)
print extract_literals(doc)
