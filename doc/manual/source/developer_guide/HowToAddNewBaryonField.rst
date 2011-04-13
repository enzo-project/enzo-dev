How to add a new baryon field
=============================

If you wish to add a new BaryonField array -- for instance, to
track the Oxygen fraction -- you'll need to do a few things.


#. Add it to the ``field_type`` structure
#. Define it in your problem type.
   :ref:`AddingANewTestProblem`.
#. Do something with it.
   :ref:`NewLocalOperator`.
#. If you need to advect a species as a color field, you will have
   to investigate how that works. Specifically, the means of
   conservation -- by fraction of by density -- as well as the inputs
   into the hydro solver.


