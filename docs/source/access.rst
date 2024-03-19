Data Access
=====

We host our public data on the Flatiron Institute computer, Rusty. Here you can find all simulations and their snapshots at redshifts 2, 1, 0.5, 0.25, and 0. If you want access to other redshifts for simulations, please get in touch with us!

.. _rusty:
CAMELS-zoomGZ on Rusty
------------
Users with a Flatiron Rusty computing account can find the hydrodynamical simulation data at:

:code:`/mnt/ceph/users/camels/PUBLIC_RELEASE/IllustrisTNG/zoom/zoom/GZ28`

This directory also hosts the file listing each simulation's astrophysical and cosmological parameters, whether they are surrogate or base simulations, and the random seed used to generate the simulations. We further host a file detailing the parameter names, fiducial values, minimum and maximum values, and a brief description.

The dark matter only analogs for the hydro simulations can be found at

:code:`/mnt/ceph/users/camels/PUBLIC_RELEASE/IllustrisTNG_DM/zoom/zoom/GZ28`


.. _binder:
CAMELS-zoomGZ on Binder
----------------
Users without a Flatiron Rusty computing account can work with the data through Binder, which is a system that allows access to the data for analysis and manipulation. More details on using Binder can be found here: `Flatiron Binder <https://wiki.flatironinstitute.org/Public/UsingFiBinder>`_. 

The binder link can be found here: `Link to Binder <https://binder.flatironinstitute.org/v2/user/sgenel/CAMELS_PUBLIC>`_.

.. warning::
    Two important things need to be taken into account when using Binder. First, the Binder environment is ephemeral - after a few days of inactivity, its contents are deleted, so one has to be vigilant
    About downloading any analysis results in time. Second, Binder is not designed to carry out long and heavy calculations. In this case, we recommend that the user download the data and work with it locally.



