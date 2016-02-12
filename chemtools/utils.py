


def doc_inherit(base_class):
    """Docstring inheriting method descriptor

       doc_inherit decorator

       Usage:

       .. code-block:: python

            class Foo(object):
                def foo(self):
                    "Frobber"
                    pass

            class Bar(Foo):
                @doc_inherit(Foo)
                def foo(self):
                    pass

       Now, ``Bar.foo.__doc__ == Bar().foo.__doc__ == Foo.foo.__doc__ ==
       "Frobber"``
    """
    def decorator(method):
        overridden = getattr(base_class, method.__name__, None)
        if overridden is None:
            raise NameError('Can\'t find method \'%s\' in base class.')
        method.__doc__ = overridden.__doc__
        return method
    return decorator
