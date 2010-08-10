Introduction to Enzo Modification
=================================

If this is the first time you've opened the hood to Enzo, welcome.
If you're an old hand and have already added new physics to it,
welcome back.

Enzo is an extremely powerful piece of software, but by no means a
complete representation of the observable universe. It's quite
likely that there will be some piece of physics that you'll want to
model, and these span a broad range of software complexities. In
all cases, whether it's a mildly invasive change such as a new
background heating model or extremely invasive like adding
relativistic non-neutral multi-fluid plasma physics, we strongly
recommend taking advantage of some basic tools. These are outlined
in the sections that follow, but the two most important ones are
`regression with LCATest? </wiki/Tutorials/LCATestSetup>`_ and
`version control with SVN? </wiki/Tutorials/MergingBranches>`_.
These two tools prevent the developer from breaking existing
features (which is far easier than one would expect), keeping track
of your changes, and sharing those changes with others. We strongly
recommend you start with getting LCATest running before you start
programming, so mistakes can be caught early.

Additionally in the Tutorials section you'll see a pair of flow
chart tools that are intended as educational tools, and several
descriptions on how to actuall add various components to the code.
It is intended that these will be at least read in order, as doing
complex things with the code require the ability to do the simpler
things.


