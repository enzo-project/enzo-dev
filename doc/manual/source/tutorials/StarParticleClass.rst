Star Particle Class
===================

Purpose
-------

To give star particles more functionality and interaction with the
grids, it was useful to create a new class for a generic particle
type
that can represent, e.g., stars, black holes, sink particles.

Main features
-------------


-  merging
-  accretion
-  conversion to a radiation source
-  adding feedback spheres to the grid, e.g. mass removal from
   accretion, supernovae.
-  different behaviors for different star types
-  multiple types of star particles
-  "active" and "inactive" stars

Approach
--------

A flowchart of the logic of the star particle class.
`Image(StarParticleClassFlowchart.pdf, 25%)? </wiki/Image(StarParticleClassFlowchart.pdf,%2025%)>`_

We keep the original implementation of the particles that are
stored
in the pointers, `ParticlePosition? </wiki/ParticlePosition>`_,
`ParticleVelocity? </wiki/ParticleVelocity>`_,
`ParticleMass? </wiki/ParticleMass>`_,
`ParticleNumber? </wiki/ParticleNumber>`_,
`ParticleType? </wiki/ParticleType>`_, and
`ParticleAttribute? </wiki/ParticleAttribute>`_. Star particles
are still created in the FORTRAN routines, e.g. star\_maker2. In
the
current version, the star class is a layer on top of these
particles.
Thus we must keep the particle pointers and objects synchronized
when
their quantities change.

Particles created in the FORTRAN routines that will be converted
into
a star object initially have a negative particle type. This
indicates
that the star is not "born" yet, which is also used to flag various
feedback spheres, such as mass removal from the grid. The stars are
activated, i.e. positive particle type, in
Star::`ActivateNewStar? </wiki/ActivateNewStar>`_ after
it has been checked for mergers, accretion, and feedback.

We store the star objects as a linked list in grid class. Because a
star object can affect multiple grids (over multiple processors)
when
adding feedback sphere, processors other than the one hosting the
star
particle needs to know about this star object. Currently for
convenience, we create a global list of star objects on all
processors. For not many stars (< 100k), this does not consume that
much memory. However in the future, we might have to reconsider how
star particles are communicated across processors.

Feedback spheres
~~~~~~~~~~~~~~~~

Any event can be set in
Star::`SetFeedbackFlag? </wiki/SetFeedbackFlag>`_ to add a feedback
sphere. This sphere can be of any size, and its properties are set
in
Star::`CalculateFeedbackParameters? </wiki/CalculateFeedbackParameters>`_
and grid::`AddFeedbackSphere? </wiki/AddFeedbackSphere>`_.
Because they can cover grids on multiple levels, we have to ensure
that they are all at the same time. In
Star::`FindFeedbackSphere? </wiki/FindFeedbackSphere>`_, we
check if sphere is completely contained within grids on the current
level. If true, we can safely add the sphere. If it's not
imperative
that the grids are completely synchronized, one can add the
feedback
sphere immediate after the star object is flagged for feedback.

Accretion / Mass Loss
~~~~~~~~~~~~~~~~~~~~~

Star objects can store up to 100 (#define MAX\_ACCR) accreti


