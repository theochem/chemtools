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
"""Test chemtools.utils.utils."""

import os
from numpy.testing import assert_equal, assert_raises
from chemtools.utils.utils import doc_inherit, Context, context


def test_doc_inherit():
    class Foo(object):
        """Dummy class for testing doc inheritance."""

        def foo(self):
            """Frobber."""
            pass

        def boo(self):
            """Boo method."""
            pass

    class Bar(Foo):
        """Dummy class for testing doc inheritance."""

        @doc_inherit(Foo)
        def foo(self):
            pass

        @doc_inherit(Foo)
        def boo(self):
            """Boo method of Bar class."""
            pass

    class Poo(Foo):
        """Dummy class for testing doc inheritance."""

        def poo(self):
            """Poo method of Poo class."""
            pass

    assert_equal(Bar.foo.__doc__, Bar().foo.__doc__)
    assert_equal(Bar.foo.__doc__, Foo.foo.__doc__)
    assert_equal(Bar.boo.__doc__, Foo.boo.__doc__)
    assert_raises(AttributeError, doc_inherit(Foo), Poo.poo)


def test_context():
    assert_equal(context.get_include(), None)
    assert_equal(hasattr(context.glob('test/*.fchk'), '__iter__'), True)
    assert_equal(context.glob('test/*.gib'), [])
    # change CTDATA
    os.environ['CTDATA'] = 'gibberish'
    assert_raises(IOError, Context)
    # remove CTDATA
    del os.environ['CTDATA']
    assert_equal(Context().data_dir[-4:], 'data')
    assert_equal(Context().get_include(), None)
