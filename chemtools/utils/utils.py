# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2014-2015 The ChemTools Development Team
#
# This file is part of ChemTools.
#
# ChemTools is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# ChemTools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
"""The Utility Module."""

import os
import sys
from glob import glob

__all__ = ['doc_inherit', 'Context', 'context']


def doc_inherit(base_class):
    """
    Docstring inheriting method descriptor.

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
        """Overwrite method docstring."""
        overridden = getattr(base_class, method.__name__, None)
        if overridden is None:
            raise NameError('Can\'t find method \'%s\' in base class.')
        method.__doc__ = overridden.__doc__
        return method

    return decorator


class Context(object):
    """
    Find out where the data directory is located etc.

    The data directory contains data files.
    """

    def __init__(self):
        # Determine data directory (also for in-place build)
        self.data_dir = os.getenv('CTDATA')
        if self.data_dir is None:
            fn_data_dir = os.path.join(os.path.dirname(__file__), 'data_dir.txt')
            if os.path.isfile(fn_data_dir):
                with open(fn_data_dir) as f:
                    self.data_dir = os.path.join(f.read().strip(), 'share/chemtools')
        if self.data_dir is None:
            self.data_dir = './data'
        self.data_dir = os.path.abspath(self.data_dir)
        # Determine include directory
        self.include_dir = os.getenv('CTINCLUDE')
        if self.include_dir is None:
            fn_data_dir = os.path.join(os.path.dirname(__file__), 'data_dir.txt')
            if os.path.isfile(fn_data_dir):
                with open(fn_data_dir) as f:
                    self.include_dir = os.path.join(
                        f.read().strip(),
                        'include/python%i.%i' % (sys.version_info.major, sys.version_info.minor))
        if not os.path.isdir(self.data_dir):
            raise IOError(
                'Can not find data files. The directory {0} does not exist.'.format(self.data_dir))

    def get_fn(self, filename):
        """Return the full path to the given filename in the data directory."""
        return os.path.join(self.data_dir, filename)

    def glob(self, pattern):
        """Return all files in the data directory that match the given pattern."""
        return glob(self.get_fn(pattern))

    def get_include(self):
        """Return the list with directories containing header files (.h and .pxd)."""
        return self.include_dir


context = Context()
