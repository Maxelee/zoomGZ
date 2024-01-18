Data Access
=====

We host our public data on the Flatiron Institute computer, Rusty. Here you can find all simulations, and their snapshots at redshifts 2, 1, 0.5, 0.25, and 0. If you want access to other redshifts for simulations please contact us!

.. _rusty:
zoomGZ on Rusty
------------
Users with a Flatiron Rusty computing account can find the hydrodynamical simulation data at:

:code:`/mnt/ceph/users/camels/PUBLIC_RELEASE/IllustrisTNG/zoom/zoomGZ`

This directory also hosts the file listing each simulation's astrophysical and cosmological parameters, whether they are surrogate or base simulations, and the random seed used to generate the simulations. We further host a file detailing the parameter names, fiducial values, minimum and maximum values, and a brief description.

The dark matter only analogs for the hydro simulations can be found at

:code:`/mnt/ceph/users/camels/PUBLIC_RELEASE/IllustrisTNG_DM/zoom/zoomGZ`


.. _binder:
zoomGZ on Binder
----------------
Users without a Flatiron Rusty computing account can work with the data through Binder, which is a system that allows access to the data for analysis and manipulation. More details on using Binder can be found here: `Flatiron Binder <https://wiki.flatironinstitute.org/Public/UsingFiBinder>`_. 

The binder link can be found here: `Link to Binder <https://binder.flatironinstitute.org/>`_.

.. warning::
    Two important things need to be taken into account when using Binder. First, the Binder environment is ephemeral - after a few days of inactivity its contents are deleted, so one has to be vigilant
    about downloading any analysis results in time. Second, Binder is not designed to carry out long and heavy calculations. In this case we recommend the user to download the data and work with it locally.



