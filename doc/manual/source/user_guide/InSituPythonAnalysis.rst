.. _in_situ_python_analysis:

In Situ Python Analysis
=======================

This feature enables Enzo to do in situ analysis using Python, by calling and using Python during simulation runtime.
This avoids outputting data to disk before conducting analysis.
Here, we use `yt <https://yt-project.org>`__ as an analysis platform.

Requirements
------------

Below are links to the build and runtime requirements, which must be installed.

* `libyt`_: a C++ shared library for in situ analysis.

  * **Normal Modes**: Shut down and terminate all the processes including simulation, if there are errors during in situ analysis using Python. This includes calling not defined functions.

  * **Interactive Modes**: Fault-tolerant to Python and supports interactive Python prompt.

* **Python >= 3.6**

  * `yt`_: an open-source, permissively-licensed python package for analyzing and visualizing volumetric data.

  * `yt_libyt`_: libyt's yt frontend.

.. _libyt: https://yt-project.github.io/libyt/HowToInstall.html#libyt

.. _yt: https://yt-project.org

.. _yt_libyt: https://yt-project.github.io/libyt/HowToInstall.html#yt_libyt

How it Works
------------
Enzo follows ``libyt``'s procedure and APIs to implement this in situ analysis feature.
You can find how ``libyt`` works more in detail `here <https://yt-project.github.io/libyt/HowItWorks.html#how-it-works>`__ and what are ``libyt`` APIs `here <https://yt-project.github.io/libyt/libytAPI>`__.

At initialization stage, ``libyt`` imports inline Python script ``inline.py`` and initializes Python interpreters in each MPI process. This happens in ``InitializeLibytInterface`` function in ``src/enzo/InitializeLibytInterface.C``.

When Enzo finishes its computation in a cycle, the whole simulation pauses and starts the in situ analysis process.
The ``CallInSitulibyt`` function in ``src/enzo/CallInSitulibyt.C`` conducts this process.
Enzo then passes in simulation information and actual field data pointers to ``libyt``.
This includes simulation information, like adaptive mesh grid hierarchy, parameters, field labels, etc, and actual simulation data pointers inside ``BaryonField`` array.  (Data fields are referenced, not copied, from ``libyt``'s perspective, although typically ``yt`` itself will make a copy as needed before any in-place changes occur.)
``libyt`` will construct data structures to store simulation information and wrap these data pointers, so that they can be read and used in Python with minimum memory overhead.
After in situ analysis is done, ``libyt`` frees resources allocated for itself, and the simulation will continue.

At the end, when the simulation is shutting down, Enzo calls ``FinalizeLibytInterface`` in ``src/enzo/InitializeLibytInterface.C``, so that ``libyt`` finalizes the Python interpreters.

How to Configure
----------------
**Some settings are hard-coded inside Enzo, you can customize it to your own needs.**

* **How to change import Python file name?**

The default Python script will be imported is ``inline.py``.
If you really want to change the name, you can go to
``src/enzo/InitializeLibytInterface.C`` in function ``InitializeLibytByItself``, and change ``params->script`` to your Python file name without ``.py``. For example, I want to make it to ``test.py``:

::

    params->script = "test";

Please use ``const char*``, or else, you have to make sure the lifetime of this variable covers the whole in situ process.

* **How to activate in situ Python analysis process?**

The full process is encapsulated inside ``CallInSitulibyt`` function in ``src/enzo/CallInSitulibyt.C``.
You can put this function everywhere you want in Enzo to start in situ analysis.
It will load and use Enzo's current state and data.

Currently, it is called inside ``EvolveLevel`` function.

* **How to call Python functions during simulation runtime? And what should I be aware of?**

You can call Python function using libyt API ``yt_run_Function`` and ``yt_run_FunctionArguments``. See how to use them `here <https://yt-project.github.io/libyt/libytAPI/PerformInlineAnalysis.html#calling-python-functions>`__.

Just put them right after the comments ``TODO: yt_run_Function and yt_run_FunctionArguments`` inside ``CallInSitulibyt`` function in ``src/enzo/CallInSitulibyt.C`` according to your needs.

Please make sure the functions you called are defined inside the script. Otherwise, in ``libyt`` normal modes, the simulation will terminate simply because it cannot find the Python function, while in interactive mode, it will labeled as failed.

See how to write an inline Python script in Doing In Situ Analysis section.

* **How to activate interactive mode and Python prompt in Enzo?**

You have to compile ``libyt`` in interactive mode.

If error occurs while running Python functions or Enzo detects ``LIBYT_STOP`` file, then ``libyt``'s interactive Python prompt will activate.

You can find more about libyt API ``yt_run_InteractiveMode`` `here <https://yt-project.github.io/libyt/libytAPI/ActivateInteractiveMode.html#activate-interactive-mode>`__.


How to Compile
--------------
The configure option that controls whether or not to use ``libyt``
can be toggled with

::

    make libyt-yes

or to turn it off,

::

    make libyt-no

1. Must use ``use-mpi-yes`` when using ``libyt-yes``. (There will be a future update to use ``libyt`` in serial, but for now, we must use with MPI.)
2. *DO NOT* use ``libyt-yes`` option and ``python-yes`` at the same time to avoid any conflicts. They are different settings.

The option will look for the following variables in the machine-specific Makefile:

::

    MACH_INCLUDES_LIBYT
    MACH_LIBS_LIBYT

If you installed ``libyt`` at ``$(LOCAL_LIBYT_INSTALL)``, which this folder include subfolders ``include`` and ``lib``, set the above variables to:

::

    MACH_INCLUDES_LIBYT = -I$(LOCAL_LIBYT_INSTALL)/include
    MACH_LIBS_LIBYT = -L$(LOCAL_LIBYT_INSTALL)/lib -lyt -Wl,-rpath,$(LOCAL_LIBYT_INSTALL)/lib

This includes ``libyt`` header, links to the library, and adds library search path for ``libyt`` library for Enzo executable.

How to Run Enzo
---------------
Put inline Python script (default file name is ``inline.py``) and Enzo executable in the same folder and run Enzo. Let's it!

If you happen to have error messages related to MPI remote memory access operation, something look like:

::

    ompi_osc_ucx_win_attach: Assertion ......

Please add ``OMPI_MCA_osc=sm,pt2pt`` before ``mpirun``, for example:

::

    OMPI_MCA_osc=sm,pt2pt mpirun -np 4 ./enzo.exe -d CollapseTestNonCosmological.enzo

This is something ``libyt`` will update and improve in the future.


Doing In Situ Analysis
----------------------
See how to write inline Python script and do in situ analysis `here <https://yt-project.github.io/libyt/InSituPythonAnalysis#in-situ-python-analysis>`__.

