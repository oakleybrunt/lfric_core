##############################################################################
# Copyright (c) 2022,  Met Office, on behalf of HMSO and the King's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################


'''PSyclone transformation script for the Dynamo0p3 API to apply the
minumum set of transformations required to permit the application to
run safely in debug mode on all supported platforms.

To reduce MPI communications, current PSyclone-LFRic strategy does not
apply halo swaps on input arguments to kernels with increment
operations on continuous fields such as GH_INC. For such kernels, PSy
layer code needs to loop into the halo to correctly compute owned dofs
on the boundary between the halo and the domain. Therefore values of
the remaining dofs in the first halo cell need to be initialised to
values that will not induce numerical errors.

By default the setval_c built-in does not initialise into
the halos. This transform causes it to do so, and so permits
developers to set safe values in halos.

'''
from psyclone.transformations import Dynamo0p3RedundantComputationTrans

def trans(psy):
    '''Applies PSyclone redundant computation transformations.

    '''
    rtrans = Dynamo0p3RedundantComputationTrans()

    setval_count = 0
    # Loop over all of the Invokes in the PSy object
    for invoke in psy.invokes.invoke_list:

        print("Transforming invoke '{0}' ...".format(invoke.name))
        schedule = invoke.schedule

        # Make setval_* compute redundantly to the level 1 halo if it
        # is in its own loop.
        for loop in schedule.loops():
            if loop.iteration_space == "dof":
                if len(loop.kernels()) != 1:
                    raise Exception(
                        "Expecting loop to contain 1 call but found '{0}'".
                        format(len(loop.kernels())))
                if loop.kernels()[0].name in ["setval_c"]:
                    setval_count += 1
                    rtrans.apply(loop, options={"depth": 1})

        # Take a look at what we've done
        print("Found {0} setval calls".format(setval_count))
        schedule.view()

    return psy
