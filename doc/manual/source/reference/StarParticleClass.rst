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

Star objects can store up to 100 (#define MAX\_ACCR) accretion
rates as
a function of time. Alternatively, currently in the black hole
particles, they can have an instantaneous accretion rate. This is
done in
Star::`CalculateMassAccretion? </wiki/CalculateMassAccretion>`_.
The actual accretion to the
star object is done in Star::Accrete.

How to add a new particle type
------------------------------

(copied from an email)


#. Set the particle type to the negative of the particle type in
   the

star maker routine. Be sure not to overwrite the type like what's
done
in the regular star\_maker routines.


2. Add the particle type to the if-statement in
   grid::`FindNewStarParticles? </wiki/FindNewStarParticles>`_.


3. Then the particles merge if any exist within

`StarClusterCombineRadius? </wiki/StarClusterCombineRadius>`_. I
should really restrict this to only star
cluster (radiating) particles. Even if there is any merging, the
particle shouldn't disappear.


4. At the end of
   `StarParticleInitialize? </wiki/StarParticleInitialize>`_, the
   routine checks if any stars

should be activated in Star\_SetFeedbackFlag. This is where I would
check first for errors or omissions. You'll have to add a new case
to
the switch statement. Something as simple as

::

    case NEW_PARTICLE_TYPE:
      if (this->type < 0)
         this->FeedbackFlag = FORMATION;
      else
         this->FeedbackFlag = NO_FEEDBACK;

will work.

After this, the particle is still negative but will be flipped
after the
feedback to the grid is applied in Star\_ActivateNewStar that's
called
from `StarParticleFinalize? </wiki/StarParticleFinalize>`_. Here
for Pop II and III stars, I use a mass
criterion. For Pop III stars, I set the mass to zero in the
pop3\_maker() f77 routine, then only set the mass after I've
applied the
feedback sphere. Perhaps you could use a similar approach... or
something more clever :)


5. The grid feedback is added in
   `StarParticleAddFeedback? </wiki/StarParticleAddFeedback>`_ that
   is

called in `StarParticleFinalize? </wiki/StarParticleFinalize>`_. In
Star\_CalculateFeedbackParameters,
you'll want to add an extra case to the switch statement that
specifies the radius of the feedback sphere and its color (metal)
density.


6. If the feedback sphere is covered by grids on the level calling

`StarParticleAddFeedback? </wiki/StarParticleAddFeedback>`_ (i.e.
all of the cells will be at the same
time), then Grid\_AddFeedbackSphere will be called. Here you'll
have to
add another if-block to add your color field to the grid.


