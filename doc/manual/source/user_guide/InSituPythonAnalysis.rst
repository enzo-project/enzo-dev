.. _in_situ_python_analysis:

In Situ Python Analysis
=======================

This feature enables Enzo to do in situ analysis using Python, by calling and using Python during simulation runtime.
This avoids outputting data to disk before conducting analysis.
Here, we use `yt <https://yt-project.org>`__ as an analysis platform.

Requirements
------------

Below are links to the build and runtime requirements, which must be installed.

``libyt`` GitHub repo can be found `here <https://github.com/yt-project/libyt>`_.
We can compile ``libyt`` using different options based on our used cases, so that Enzo can have different in situ analysis feature when it links to ``libyt``.

A brief description of each mode (option) is shown here. The options are for compiling ``libyt`` only, and they are mutually independent.
Please follow the instructions in ``libyt`` documents:

* `libyt`_: a C shared library for in situ analysis.

  * **Serial Mode** (``-DSERIAL_MODE=ON``): Compile ``libyt`` using GCC compiler.

  * **Parallel Mode** (``-DSERIAL_MODE=OFF``): Compile ``libyt`` using MPI compiler and use it in parallel computation.

  * **Normal Mode** (``-DINTERACTIVE_MODE=OFF``): Shut down and terminate all the processes including simulation, if there are errors during in situ analysis using Python.

  * **Interactive Mode** (``-DINTERACTIVE_MODE=ON``): Fault-tolerant to Python and support interactive Python prompt and reloading script feature.

  * **Jupyter Kernel Mode** (``-DJUPYTER_KERNEL=ON``): Fault-tolerant to Python and support using Jupyter Notebook / JupyterLab UI to connect to simulation and do in situ analysis.

* A Python >= 3.7 environment with the following packages installed:

  * `mpi4py`_: It provides Python bindings for the Message Passing Interface (MPI) standard. It is required for **parallel mode**. ``mpi4py`` installed should match the MPI used for compiling Enzo.

  * `yt`_: An open-source, permissively-licensed python package for analyzing and visualizing volumetric data.

  * `yt_libyt`_ (>=0.0.9): ``libyt``'s yt frontend.

  * `jupyter_libyt`_: A Jupyter Client plugin for connecting to libyt Jupyter kernel. This is only required in **jupyter kernel mode**.

.. _libyt: https://libyt.readthedocs.io/en/latest/

.. _yt: https://yt-project.org

.. _yt_libyt: https://github.com/data-exp-lab/yt_libyt

.. _jupyter_libyt: https://github.com/yt-project/jupyter_libyt

.. _mpi4py: https://mpi4py.readthedocs.io/en/stable/install.html#installation

How it Works
------------
Enzo follows ``libyt``'s procedure and APIs to implement this in situ analysis feature.
We can find how ``libyt`` works more in detail `here <https://libyt.readthedocs.io/en/latest/how-it-works.html>`__ and what are ``libyt`` APIs `here <https://libyt.readthedocs.io/en/latest/libyt-api/index.html>`__.

At initialization stage, ``libyt`` imports a Python script (the default name is ``inline.py``; see below for using alternative names) and initializes Python interpreters in each MPI process (if compiled with ``-DSERIAL_MODE=OFF``). This happens in ``InitializeLibytInterface`` function in ``src/enzo/InitializeLibytInterface.C``.

When Enzo finishes its computation in a cycle, the whole simulation pauses and starts the in situ analysis process.
The ``CallInSitulibyt`` function in ``src/enzo/CallInSitulibyt.C`` conducts this process.
Enzo then passes in simulation information and actual field data pointers to ``libyt``.
This includes simulation information, like adaptive mesh grid hierarchy, parameters, field labels, etc, and actual simulation data pointers inside ``BaryonField`` array.  (Data fields are referenced, not copied, from ``libyt``'s perspective, although typically ``yt`` itself will make a copy as needed before any in-place changes occur.)
``libyt`` will construct data structures to store simulation information and wrap these data pointers, so that they can be read and used in Python with minimum memory overhead.

``libyt`` supports calling Python functions from simulation process,
and it also supports user interface (**interactive mode** has :ref:`Interactive Python Prompt` and :ref:`Reloading Script` feature, and **jupyter kernel mode** has :ref:`Jupyter Notebook / JupyterLab UI` feature).
We can update Python functions and probe data interactively in the UI.
After in situ analysis is done, ``libyt`` frees resources allocated for itself, and the simulation will continue.

At the end, when the simulation is shutting down, Enzo calls ``FinalizeLibytInterface`` in ``src/enzo/InitializeLibytInterface.C``, so that ``libyt`` finalizes the Python interpreters.

How to Configure
----------------
**Some settings are hard-coded inside Enzo, we can customize it to our own needs.**

General
^^^^^^^

* **How to change import Python file name?**

  The default Python script will be imported is ``inline.py``.

  If we want to change the script name, set the file name without its extension to ``libyt_script_name`` in Enzo parameter file.
  For example, we want to make the Python script to be ``test.py``:

  ::

      libyt_script_name = test


* **How to call Python functions during simulation runtime? And what should I be aware of?**

  We can call Python function using libyt API ``yt_run_Function`` and ``yt_run_FunctionArguments``. See how to use them `here <https://libyt.readthedocs.io/en/latest/libyt-api/run-python-function.html>`__.

  Un-comment the code in ``src/enzo/CallInSitulibyt.C``:

  ::

    /* Run yt_run_Function and yt_run_FunctionArguments */
    // if (yt_run_Function("yt_inline") != YT_SUCCESS) {
    // 	   fprintf(stderr, "Error while running yt_run_Function and call yt_inline\n");
    // 	   return FAIL;
    // }
    // if (yt_run_FunctionArguments("yt_inline_args", 1, "\'density\'") != YT_SUCCESS) {
    //     fprintf(stderr, "Error while running yt_run_FunctionArguments and call yt_inline_args\n");
    //     return FAIL;
    // }

  The corresponding inline script ``inline.py`` would be:

  ::

    import yt_libyt
    import yt
    #yt.enable_parallelism() # make yt works in parallel computing (require mpi4py)

    def yt_inline():
        pass

    def yt_inline_args(field):
        pass

  Please make sure the functions we called are defined inside the script. Otherwise, in ``libyt`` normal modes, the simulation will terminate simply because it cannot find the Python function, while in the other modes, it will labeled as failed.

  See how to use yt to do analysis `here <https://libyt.readthedocs.io/en/latest/in-situ-python-analysis/using-yt.html>`__.

.. _Interactive Python Prompt:

Interactive Python Prompt
^^^^^^^^^^^^^^^^^^^^^^^^^

* **How to activate interactive Python prompt in Enzo?**

  We have to compile ``libyt`` in **interactive mode** and then use ``libyt-yes`` and ``libyt-interactive-yes`` options to compile Enzo.

  The code in ``src/enzo/CallInSitulibyt.C`` will call libyt API:

  ::

    /* Call interactive Python prompt. */
    if (yt_run_InteractiveMode("LIBYT_STOP") != YT_SUCCESS) {
        fprintf(stderr, "Error in libyt API yt_run_InteractiveMode\n");
        fprintf(stderr, "One reason might be compiling libyt without -DINTERACTIVE_MODE=ON, "
                        "which does not support yt_run_InteractiveMode.\n");
    }

  If Enzo detects ``LIBYT_STOP`` file, then interactive Python prompt will activate.
  We can find more about libyt API ``yt_run_InteractiveMode`` `here <https://libyt.readthedocs.io/en/latest/libyt-api/yt_run_interactivemode.html>`__.

* **How to use interactive Python prompt? How does it work?**

  It is like a normal Python prompt but with access to simulation data,
  see `here <https://libyt.readthedocs.io/en/latest/in-situ-python-analysis/interactive-python-prompt.html>`__ for how to use interactive Python prompt.

  Interactive Python prompt only works on local desktops or submit an interactive job to HPC cluster (ex: ``qsub -I`` in PBS scheduler),
  because the prompt gets inputs from the terminal.
  The root process gets the inputs and then broadcasts the inputs to every MPI process. They run the statements synchronously.

.. _Reloading Script:

Reloading Script
^^^^^^^^^^^^^^^^

* **How to activate reload Python script in Enzo?**

  We have to compile ``libyt`` in **interactive mode** and then use ``libyt-yes`` and ``libyt-reload-yes`` options to compile Enzo.

  The code in ``src/enzo/CallInSitulibyt.C`` will call libyt API:

  ::

    /* Reloading script */
    if (yt_run_ReloadScript("LIBYT_STOP", "RELOAD", "reload.py") != YT_SUCCESS) {
        fprintf(stderr, "Error in libyt API yt_run_ReloadScript\n");
        fprintf(stderr, "One reason might be compiling libyt without -DINTERACTIVE_MODE=ON, "
                        "which does not support yt_run_ReloadScript.\n");
    }

  If an error occurred when running inline Python functions or Enzo detects ``LIBYT_STOP`` file, then it will enter reloading script stage.
  Document about ``yt_run_ReloadScript`` is `here <https://libyt.readthedocs.io/en/latest/libyt-api/yt_run_reloadscript.html>`__.

* **How to reload a script?**

  Reloading script feature is a file-based interactive Python prompt, such that user creates specific files to send instructions to libyt and gets outputs from a file.
  The feature can be used in HPC cluster and does not limit to interactive jobs only.

  Some used cases are, for example, when an unexpected Python error occurred during the simulation runtime, we can update the function just in time and do not need to go all over again;
  or when we want to change the Python script during runtime.

  See `here <https://libyt.readthedocs.io/en/latest/in-situ-python-analysis/reloading-script.html>`__ for how to reload a script.

.. _Jupyter Notebook / JupyterLab UI:

Jupyter Notebook / JupyterLab UI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **How does this work?**

  Traditionally, when we start a Jupyter Notebook, it summons a kernel and manages it itself.
  A kernel is the actual core process that runs the commands and statements in a Jupyter Notebook (JupyterLab) cell.

  Here, we do it in an opposite way.
  We launch a kernel first, and then we make Jupyter Notebook connects to it. The kernel is no longer managed by Jupyter at all.

  We run simulation and start Jupyter Notebook separately.
  libyt API ``yt_run_JupyterKernel`` launches a libyt Jupyter kernel (libyt kernel), so that simulation data is within reach.
  We then use Jupyter Notebook to connect to libyt kernel.
  Thus we can do in situ analysis using Jupyter Notebook UI.

* **How to launch libyt kernel in Enzo?**

  We have to compile ``libyt`` in **jupyter kernel mode** and then use ``libyt-yes`` and ``libyt-jupyter-yes`` options to compile Enzo.

  The code in ``src/enzo/CallInSitulibyt.C`` will call libyt API:

  ::

    /* Launch libyt Jupyter kernel */
    if (yt_run_JupyterKernel("LIBYT_STOP", false) != YT_SUCCESS) {
         fprintf(stderr, "Error in libyt API yt_run_JupyterKernel\n");
         fprintf(stderr, "One reason might be compiling libyt without -DJUPYTER_KERNEL=ON, "
                         "which does not support yt_run_JupyterKernel.\n");
    }

  If Enzo detects ``LIBYT_STOP`` file, it will launch a libyt kernel.
  Since we set the second argument to ``false``, libyt kernel will bind to empty ports automatically.
  If it is set to ``true``, libyt kernel will use the configuration based on user-provided connection file.
  This is useful when running simulation in HPC clusters.
  See `here <https://libyt.readthedocs.io/en/latest/libyt-api/yt_run_jupyterkernel.html>`__ for libyt API ``yt_run_JupyterKernel``.

* **How to start Jupyter Notebook / JupyterLab and connect to libyt kernel? How to use it?**

  This feature can be used in local desktop and HPC cluster.
  See `here <https://libyt.readthedocs.io/en/latest/in-situ-python-analysis/jupyter-notebook/jupyter-notebook-access.html>`__ for a step by step guide and how to use it.

  Notice that ``libyt`` hasn't done implementing Jupyter's full feature.
  What it does is processing inputs and printing outputs faithfully.
  Features like data streaming, debugging, and ipwidgets are not supported yet.
  ``libyt`` will add these features in the future update.

How to Compile
--------------
The configure option that controls whether or not to use ``libyt`` can be toggled with:

::

    make libyt-yes

or to turn it off,

::

    make libyt-no

There are other subsettings, and they must use ``libyt-yes``:

* ``libyt-interactive-yes``/``libyt-interactive-no``: set whether to use libyt :ref:`Interactive Python Prompt`

* ``libyt-reload-yes``/``libyt-reload-no``: set whether to use libyt :ref:`Reloading Script` feature

* ``libyt-jupyter-yes``/``libyt-jupyter-no``: set whether to use :ref:`Jupyter Notebook / JupyterLab UI`.

*DO NOT* use ``libyt-yes`` option and ``python-yes`` at the same time to avoid any conflicts. They are different settings.

The option ``libyt-yes`` will look for the following variables in the machine-specific Makefile:

::

    MACH_INCLUDES_LIBYT
    MACH_LIBS_LIBYT

If we installed ``libyt`` at ``$(LOCAL_LIBYT_INSTALL)``, which this folder includes subfolders ``include`` and ``lib``, set the above variables to:

::

    MACH_INCLUDES_LIBYT = -I$(LOCAL_LIBYT_INSTALL)/include
    MACH_LIBS_LIBYT = -L$(LOCAL_LIBYT_INSTALL)/lib -lyt -Wl,-rpath,$(LOCAL_LIBYT_INSTALL)/lib

This includes ``libyt`` header, links to the library, and adds library search path for ``libyt`` library for Enzo executable.

How to Run Enzo
---------------
Put inline Python script (default file name is ``inline.py``) and Enzo executable in the same folder and run Enzo.

If we happen to have error messages related to MPI remote memory access operation, something look like:

::

    ompi_osc_ucx_win_attach: Assertion ......

Please add ``OMPI_MCA_osc=sm,pt2pt`` before ``mpirun``, for example:

::

    OMPI_MCA_osc=sm,pt2pt mpirun -np 4 ./enzo.exe -d CollapseTestNonCosmological.enzo

It is for one-sided MPI communication settings.

.. _Doing In Situ Analysis:

Doing In Situ Analysis
----------------------

Generally, after Enzo has loaded simulation data to Python, ``libyt`` can analyze data using arbitrary Python script in parallel computing using ``mpi4py``.
Every Python instance synchronously runs the statement line-by-line inside the imported script's namespace.
If there are *N* simulation processes, then there will be *N* Python instances working together to conduct in situ analysis.

``yt`` already supports `parallel computation <https://yt-project.org/doc/analyzing/parallel_computation.html#parallel-computation-with-yt>`__,
so we can directly use it in the Python script.
When converting post-processing script to inline Python script, remember to:

* import ``yt_libyt`` and change ``yt.load`` to ``yt_libyt.libytDataset()``
* call ``yt.enable_parallelism()`` if Enzo is running in MPI platform.

For example, this creates a radial profile of the density field at the domain center:

::

  import yt_libyt
  import yt
  #yt.enable_parallelism()   # make yt works in parallel computing in libyt parallel mode (require mpi4py)

  def yt_inline():
      ds = yt_libyt.libytDataset()
      sphere = ds.sphere(ds.domain_center, (20.0, "km"))
      profile = yt.ProfilePlot(sphere, ("index", "radius"), ("gas", "density"))

      # save the figure on root process only
      if yt.is_root():
          profile.save()

See how to write inline Python script using ``yt`` `here <https://libyt.readthedocs.io/en/latest/in-situ-python-analysis/using-yt.html>`__.
