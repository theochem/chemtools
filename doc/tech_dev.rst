..
    : ChemTools is a collection of interpretive chemical tools for
    : analyzing outputs of the quantum chemistry calculations.
    :
    : Copyright (C) 2016-2019 The ChemTools Development Team
    :
    : This file is part of ChemTools.
    :
    : ChemTools is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : ChemTools is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --


.. _usr_development:

Developer Guidelines
####################

Submitting Issues
=================

In case you experience issues or bugs using this package, please feel free to
submit an issue on our `issue page <https://github.com/QuantumElephant/chemtools/issues>`_.

.. _dev_build:

Building ChemTools
==================

We use github to host our repository. You can download our
`latest sources <https://github.com/QuantumElephant/chemtools>`_ with command:

.. code-block:: bash

   git clone https://github.com/QuantumElephant/chemtools.git chemtools

For developer, it is more convenient to build ChemTools in place. Before you build
ChemTools, ensure you have all the required :ref:`dependencies <usr_py_depend>` installed.
You can build ChemTools with the commands:

.. code-block:: bash

   pip install -e .[dev]


Running Tests
=============

Nosetests are recommended after built to ensure all the parts of ChemTools
working properly:

.. code-block:: bash

   nosetests -v chemtools

At this stage, some ``UserWarning`` messages are printed in between tests which is expected.
However, no test should fail.

Contributing Code
=================

The preferred way to contribute to ChemTools is fork from our main repo
on github, make some changes, and then make a `Pull request`:

* Fork the `main repo <https://github.com/QuantumElephant/chemtools>`_ from ChemTools github page: Clicking on the ``Fork`` button
  near the top right of the page. This will create a copy of the latest code
  under your own account.

* Add your forked repo to the remote list of your local repo:

  .. code-block:: bash

     git remote add your_remote_name git@github.com:YourAccount/chemtools.git

* Create a branch to hold your changes:

  .. code-block:: bash

     git checkout -b your_feature_branch_name

* Include your name and email in the ``AUTHORS`` file at the root directory of ChemTools.
  Ensure the name and email are the same as your git config. You can check your
  name and email with the commands:

  .. code-block:: bash

     git config user.name
     git config user.email

* Work on this branch locally and make some changes. When you've done editing,
  use the command:

  .. code-block:: bash

     git add files_you_modified
     git commit

  to write down the changes your made and save it.

* Push your changes to your own remote repo:

  .. code-block:: bash

     git push your_remote_name your_feature_branch_name

* Finally, go to your forked github repo page, click ``Pull request`` to send your
  changes. All the changes need to pass the automatic quality test before your
  pull request gets reviewed. You can go to the `"Pull request" <https://github.com/QuantumElephant/chemtools/pulls>`_ page of the main repo
  to check the status of the test and fix the errors if any of them fail.

.. _usr_doc:

Building Documentation
======================

If you are interested in generating the documentation from source, the following
packages are also needed:

* Sphinx >=1.3.1: http://sphinx.pocoo.org/
* sphinxcontrib-bibtex >= 0.3.5: https://pypi.python.org/pypi/sphinxcontrib-bibtex
* IPython >= 3.2.1: https://ipython.org/install.html

To install these dependencies,

.. code-block:: bash

   pip install -e .[doc]



The Sphinx Read-The-Docs theme as well as Sphinx Gallery customized for ChemTools can be
obtained by cloning the repository as a submodule from ChemTools parent directory:

.. code-block:: bash

    git submodule update --init --recursive


.. Also, make sure that the :ref:`examples files are downloaded <usr_lfs_installation>`.

To automatically generate API documentation and generate HTML:

.. code-block:: bash

   cd doc
   make clean
   make html

To open the documentation in your default browser, either click on ``_build/html/index.html``
file directly, or run the command below from terminal:

.. code-block:: bash

   open _build/html/index.html

In case this command did not work, for example on Ubuntu 16.04 you may get a message like **"Couldn't get a
file descriptor referring to the console"**, try:

.. code-block:: bash

   see _build/html/index.html


Quality Assurance
=================

When making a pull request to contribute to the ChemTools repository, the code is remotely tested to see
if it passes all the tests and meets ChemTools' quality standards. To run the tests locally, please refer
to :ref:`Testing <usr_testing>`. If you are interested to run the quality assurance scripts locally, first
install the dependencies below:

* PyLint >= 1.5.0: https://www.pylint.org/
* pycodestyle >= 2.0.0: http://pycodestyle.readthedocs.io/
* pydocstyle >= 1.0.0: http://pydocstyle.readthedocs.io/
* coverage >= 4.1: https://coverage.readthedocs.io/
* Git >= 1.8: https://git-scm.com/
* GitPython >= 2.0.5: http://gitpython.readthedocs.io/

Then, download the quality assurance code by cloning the corresponding submodule:

.. code-block:: bash

   git submodule update --init tools/inspector

And, run the module's bash script to setup some pre-commit hooks and copy files to run the quality assurance
scripts individually:

.. code-block:: bash

   # it is installed in the relative path
   cd tools/inspector
   ./install.sh
   cd ../..

At this stage, the quality assurance tests can be simulated from the ChemTools parent directory.
For example to run ``pylint`` check,

.. code-block:: bash

   # from ChemTools parent directory
   ./tools/inspector/qa/simulate_trapdoor_pr.py tools/inspector/qa/trapdoor_pylint.py

To run all of the quality assurance scripts,

.. code-block:: bash

   # from ChemTools parent directory
   for i in tools/inspector/qa/trapdoor_*; do tools/inspector/qa/simulate_trapdoor_pr.py $i; done

Note that you should be developing on a feature (not master) branch and merging/rebasing to the
updated master when complete. There should be also no uncommitted changes when running these scripts.
