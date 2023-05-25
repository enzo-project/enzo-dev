.. _in_situ_python_analysis:

In Situ Python Analysis
=======================
This feature enables Enzo to do in situ analysis using Python, which is to call and use Python during simulation runtime.
We no longer need to dump data on hard disk first, before we do any analysis.
Here, we use `yt <http://yt-project.org>`_ as its core analytic method.

Requirements
------------
Please follow the links to install the requirements.

* **`libyt <https://calab-ntu.github.io/libyt/HowToInstall.html#libyt>`_**: a C++ library for in situ analysis.
* Python >= 3.6
  * **`yt <http://yt-project.org>`_**: an open-source, permissively-licensed python package for analyzing and visualizing volumetric data.
  * **`yt_libyt <https://calab-ntu.github.io/libyt/HowToInstall.html#yt_libyt>`_**: yt frontend for libyt.

How it Works
------------
Enzo follows ``libyt``'s procedure and APIs to implement this in situ analysis feature.
You can find how ``libyt`` works more in detail `here <https://calab-ntu.github.io/libyt/HowItWorks.html#how-it-works>`_ and what are ``libyt`` APIs `here <https://calab-ntu.github.io/libyt/libytAPI>`_.

At initialization stage, ``libyt`` imports inline Python script ``inline.py`` and initializes Python interpreters in each MPI process. This happens in ``InitializeLibytInterface`` function in ``src/enzo/InitializeLibytInterface.C``.

When finishes its computation in a cycle, the whole simulation pauses and starts the in situ analysis process.
``CallInSitulibyt`` function in ``src/enzo/CallInSitulibyt.C`` goes through the whole process.
It passes in simulation information and actual field data pointers to ``libyt``.
This includes simulation information, like adaptive mesh grid hierarchy, parameters, field labels, etc, and actual simulation data pointers inside ``BaryonField`` array.
``libyt`` will construct data structure to store simulation information and wrap these data pointers, so that they can be read and used in Python with minimum memory overhead.
After in situ analysis is done, ``libyt`` frees resources allocated for itself, and the simulation will continue.

At the end, when the simulation is shutting down, Enzo calls ``FinalizeLibytInterface`` in ``src/enzo/InitializeLibytInterface.C``, so that ``libyt`` finalizes the Python interpreters.

How to Configure
----------------


How to Compile
--------------


How to Run Enzo
---------------


Doing In Situ Analysis
----------------------


Frequently Asked Questions
--------------------------
* **Problems when RMA...?**

* **What if I get segmentation fault?**
