How to add a new baryon field
=============================

If you wish to add a new BaryonField array -- for instance, to
track the Oxygen fraction -- you'll need to do a few things.


#. Add it to the ``field\_type`` structure
#. Define it in your problem type.
   `(You'll probably also add a new test problem)? </wiki/Tutorials/NewTestProblem1>`_
#. Do something with it.
   `(You'll probably need to add a new local operator, as well.)? </wiki/Tutorials/NewLocalOperator>`_
#. If you need to advect a species as a color field, you will have
   to investigate how that works. Specifically, the means of
   conservation -- by fraction of by density -- as well as the inputs
   into the hydro solver.


