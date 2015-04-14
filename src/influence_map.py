import sys
import os
import re
from pysb import kappa

def main(argv):
    if len(sys.argv) != 3:
        print "Usage: %s model_filename im_filename" % __file__
        return 1

    model_filename = argv[1]
    im_filename = argv[2]

    # Sanity checks on filename
    if not os.path.exists(model_filename):
        raise Exception("File '%s' doesn't exist" % model_filename)
    if not re.search(r'\.py$', model_filename):
        raise Exception("File '%s' is not a .py file" % model_filename)
    sys.path.insert(0, os.path.dirname(model_filename))
    model_name = re.sub(r'\.py$', '', os.path.basename(model_filename))
    # import it
    try:
        # FIXME if the model has the same name as some other "real" module
        # which we use, there will be trouble (use the imp package and import
        # as some safe name?)
        model_module = __import__(model_name)
    except StandardError as e:
        print "Error in model script:\n"
        raise
    # grab the 'model' variable from the module
    try:
        model = model_module.__dict__['model']
    except KeyError:
        raise Exception("File '%s' isn't a model file" % model_filename)

    # Make the contact map
    im_output_filename = kappa.influence_map(model, base_filename=im_filename)

    # Process the influence map file to fix the extra newlines
    im_output_file = open(im_output_filename)
    im_output_file_fixed = open('%s_im_fixed.dot' % im_filename, 'w')
    for line in im_output_file.readlines():
        fixed_line = re.sub('label="\n', 'label="', line)
        im_output_file_fixed.write(fixed_line)
    im_output_file.close()
    im_output_file_fixed.close()
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))

