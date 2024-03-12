from inspect import signature
from dtcc_model.model import Model


def register_model_method(fn):
    params = signature(fn).parameters
    if len(params) == 0:
        raise ValueError("Model method must have at least one parameter")
    first_arg = list(params.keys())[0]
    first_type = params[first_arg].annotation
    if not issubclass(first_type, Model):
        raise ValueError(
            f"First parameter must be a DTCC Model. Did you forget to add a type hint?"
        )
    fn._dtcc_model_method = True
    first_type.add_methods(fn, fn.__name__)
    return fn
