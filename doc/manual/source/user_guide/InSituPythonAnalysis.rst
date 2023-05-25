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
At initialization, Enzo imports inline Python script ``inline.py``. (START HERE)

When finishes its computation in a cycle, the whole simulation pauses and starts the in situ analysis process.
Enzo passes in simulation information and actual field data pointers to ``libyt``.
This includes adaptive mesh grid hierarchy, parameters, field labels, data pointers inside ``BaryonField`` array, etc.
``libyt`` will construct data structure to store

You can find more about how ``libyt`` works `here <https://calab-ntu.github.io/libyt/HowItWorks.html#how-it-works>`_.


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
