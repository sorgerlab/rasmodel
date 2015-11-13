from pysb.util import rules_using_parameter
from pysb import Parameter, ComponentSet

def merge_parameters(model, new_name, parameters):
    unique_values = {parameter.value for parameter in parameters}
    if len(unique_values) > 1:
        raise ValueError("Given parameters have different values: %s" %
                         (', '.join('%s=%g' % (p.name, p.value)
                                    for p in parameters)))
    value = parameters[0].value
    rules = ComponentSet()
    for parameter in parameters:
        rules |= rules_using_parameter(model, parameter)
    if not rules:
        raise ValueError("Model has no rules using given parameters: %s" %
                         ', '.join(p.name for p in parameters))
    try:
        new_parameter = model.parameters[new_name]
        if new_parameter.value != value:
            raise ValueError("Parameter %s is already present in the model "
                             "with the value %g, which differs from the "
                             "common value of the given parameters, %g" %
                             (new_parameter.name, new_parameter.value, value))
    except KeyError:
        new_parameter = Parameter(new_name, value)
        model.add_component(new_parameter)
    for rule in rules:
        for attr in 'rate_forward', 'rate_reverse':
            if getattr(rule, attr) in parameters:
                setattr(rule, attr, new_parameter)
