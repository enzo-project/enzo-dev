.. _BuildingDocumentation:

Building the Documentation
============================

The documentation for Enzo (including this very document) is built using
the `reStructuredText (ReST) <http://docutils.sourceforge.net/rst.html>`_
syntax which is parsed into final formats using the
`Sphinx engine <http://sphinx.pocoo.org/>`_.
Sphinx is a python package which may be installed using the 
`pip <http://www.pip-installer.org/en/latest/>`_ Python package installation
tool like this:

.. highlight:: none

::

  $ pip install sphinx

Once that is installed, make sure that the binary ``sphinx-build`` is in your
path (``$ which sphinx-build``). Relative to the top level of the Enzo package,
the Enzo docs are in ``doc/manual``.
This directory contains a ``Makefile`` and a ``source`` directory.
From within this directory, this command will parse the documents into
a hierarchy of HTML files (identical what is on the web) into a
new directory ``build``:

::

  $ make clean
  $ make html

If that is successful, point your web browser to the file on disk
(using the Open File... option of the File menu) ``build/html/index.html``
(this is relative to this same directory with the ``Makefile``).
On Mac OS X this command should work: ``open build/html/index.html``.
The docs should be nearly
identical to what is online, but they are coming from the local machine.


Building a PDF of the Documentation
------------------------------------

If (PDF)LaTeX is functional, is it possible to build a PDF of the Enzo
documentation in one step.
In the directory with the ``Makefile``, use this command:

::

  $ make latexpdf

If this is successful, the PDF will be ``build/latex/Enzo.pdf``.
The PDF might be preferred for some users, and can be searched all at once for
a term, unlike a local copy of the HTML.


If PDFLaTeX
is not working, ``$ make latex`` will not attempt to make the PDF. A PS or DVI
(or whatever anachronistic thing your SPARCstation makes)
can be made starting from ``build/latex/Enzo.tex``.


Updating the Online Pre-Built Documentation
---------------------------------------------

If you are an Enzo developer and need to update the current build of the
documentation, these instructions should step you through that process.
In the directory doc/manual, clone a copy of the
`docs repository <http://code.google.com/p/enzo/source/checkout?repo=docs>`_
and call it enzo-docs:

::

  $ hg clone https://docs.enzo.googlecode.com/hg/ enzo-docs

The enzo-docs repository holds HTML and image files which are accessed
over the web on the Google Code pages, and it is not
modified manually.

In the directory doc/manual, execute these commands:

::

  $ cd enzo-docs
  $ hg pull
  $ hg up
  $ cd ..
  $ rm -rf enzo-docs/*
  $ make clean
  $ make html
  $ export CURRENT_REV="`hg identify`"
  $ cp -Rv build/html/* enzo-docs/
  $ cd enzo-docs
  $ hg addremove --similarity=25
  $ hg ci -m "Build from ${CURRENT_REV}"
  $ hg push

The ``addremove`` extension will automatically update the current state of the
files to be pushed with whatever is found in that directory. It will remove all
files that it no longer sees and add all the files it currently does see.
(You can supply ``--dry-run`` to see what it will do). The similarity argument just
helps with keeping the size of the commits down, but because this repository
is only to be used as a holding place for static content this should not be a
problem.

