.. _DoingARelease:

Doing a Release
===============

Periodically, the Enzo community creates an "official" release of the Enzo
codebase.  While running off of the Enzo mercurial repository is in general
quite stable, doing official releases has other benefits.  By doing releases, we
acoomplish the following goals:

* Periodically recognize the breadth and depth of the code contributions by our
  community.
* Offer a "stable" platform with known properties to test against for people who
  are not heavily involved in Enzo development.
* Announce to the wider computational astrophysics community about ongoing
  developments in the Enzo codebase.

Generally, releases happen via the contributions of a release manager and the
author of the release e-mail.  

Generally, the release manager is a senior member of the community whose
responsibility is to ensure open pull requests are integrated into the code
before the release, select a release e-mail author, and ensure that the
checklist in this document is carried out.

The author of the release e-mail is generally someone who has made significant
recent contributions to the code.  This person may be at any seniority level,
although in the past several releases (as of Enzo 2.4) this person has generally
been either a postdoc or a grad student.

To do the release, the following tasks must be completed:

* Update the ``README`` file in the root of the repository to reflect the
  current version. Also look over the document to correct any changes to
  repository locations, mailing list or social media addresses, or new
  contributors.

* Update the ``CHANGELOG`` to include a new entry for the release.  The
  demarcation between new features, enhancements, or bugfixes is up to the
  judgement of the release manager. Use the following format::

   == Version 2.x ==
   _Release Date: 1/19/2038
 
   * New Feature: A frobulator was added to the code to improve frobulation.
                  (PR xxx)
   * Enhancement: The moving mesh module now supports 11-dimensional meshes.
                  (PR YYY)
   * Bugfix: The retro-encabulator no longer instantiates sentient AIs. 
             (PR ZZZ)
  
* Update the ``conf.py`` file in the documentation to include the new version
  number.

* Ensure that the answer tests are passing on the automated build machine.

* Once all pull requests slated for the release have been merged, tag the final
  commit as the "enzo-2.x" release changeset.
