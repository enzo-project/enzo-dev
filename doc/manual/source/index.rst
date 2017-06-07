.. Enzo documentation master file, created by
   sphinx-quickstart on Fri Jul  2 15:20:02 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Enzo's documentation!
================================

This is the development site for Enzo, an adaptive mesh refinement (AMR),
grid-based hybrid code (hydro + N-Body) which is designed to do simulations of
cosmological structure formation. Links to documentation and downloads for all
versions of Enzo from 1.0 on are available. 

Enzo development is supported by grants AST-0808184 and OCI-0832662 from the
National Science Foundation. 


.. toctree::
   :maxdepth: 2

   EnzoLicense.rst
   user_guide/index.rst
   supplementary_info/index.rst
   parameters/index.rst
   physics/index.rst
   faq/index.rst
   developer_guide/index.rst
   reference/index.rst
   presentations/index.rst

Enzo Mailing Lists
-----------------------

There are two mailing lists for Enzo hosted on Google Groups, enzo-users
and enzo-dev.

enzo-users
^^^^^^^^^^

Everyone Enzo user should sign up for the
enzo-users mailing list.
This is is used to announce changes to Enzo, and sometimes major changes
Enzo-related analysis tools.
This list is appropriate for anything else Enzo-related,
such as machine-specific compile problems,
discussions of the science and physics behind what Enzo does,
or queries about problem initialization.
We recommend using the Enzo users mailing list liberally - by this we mean
that any question asked on the list will educate everyone else on the list,
and is manifestly not a stupid question.
As long as a good effort has been made to try to figure out the answer
before mailing the list, all questions about Enzo are welcome!
Please follow the link below to sign up for this list and a link
to discussion archives:

`http://groups.google.com/group/enzo-users <http://groups.google.com/group/enzo-users>`_

To post a message to this list, send an email to:

enzo-users@googlegroups.com

The archives for the old Enzo users mailing list can be found linked below.
A search of the list archives should be performed before emailing the list
to prevent asking a question that has already been answered (using, for example,
`an advanced web search <http://www.google.com/advanced_search>`_
limited to that page).

https://mailman.ucsd.edu/pipermail/enzo-users-l/

enzo-dev
^^^^^^^^

The second mailing is for developers of Enzo. This is for Enzo "old-hats",
or anyone interested in adding new features to Enzo, or anyone who wants a deeper
understanding of the internals of Enzo. Please follow the link below
to sign up for the list and a link to the discussion archives:

`http://groups.google.com/group/enzo-dev <http://groups.google.com/group/enzo-dev>`_

To post a message to this list, send an email to:

enzo-dev@googlegroups.com



Regression Tests
----------------

Enzo has an internal testing suite (:ref:`EnzoTestSuite`) that
performs regression tests that verifies that the code is producing
expected results on a wide variety of platforms.  It also aids in
discovering bugs that may have been introduced in the development
process of Enzo.  The Enzo codebase is tested before every point
release and routinely by Enzo developers.

Citing Enzo
-----------

Guidelines for citing enzo are available in the ``CITATION`` file in the root of
the enzo mercurial repository.

If you use Enzo for a scientific publication, we ask that you cite the code in
the following way in the acknowledgments of your paper::

    Computations described in this work were performed using the
    publicly-available \texttt{Enzo} code (http://enzo-project.org), which is
    the product of a collaborative effort of many independent scientists from
    numerous institutions around the world.  Their commitment to open science
    has helped make this work possible.

In addition, we request that you link to the project webpage in a footnote and
add a citation to the Enzo method paper.  See the ``CITATION`` file for BibTeX
and LaTeX formatted citations.

Search
======

* :ref:`search`

