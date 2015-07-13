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
documentation, simply modify the docs in the enzo-dev repository in the
same way you would edit the source code.  The docs exist in the 
enzo-dev/doc directory.  Submit a pull request for these changes in the
same way you would do so with source modifications.  If accepted, 
these new docs will be available almost immediately at:
``http://enzo.readthedocs.org``.
