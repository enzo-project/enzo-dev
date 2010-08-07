Using Parallel Root Grid IO
===========================

First, read
`How Does Parallel Root Grid IO Work.? </wiki/HowDoesParallelRootGridIOwork>`_
Come back when you're finished.

Your PRGIO problem generator needs to function in two passes. It
needs to set up the basic problem (see
`Part 1: Serial Unigrid? </wiki/NewTestProblem/Part1_SerialUnigrid>`_)
*without* allocating any data. Then, after
``CommunicationPartitionGrid`` is called, it needs to *re*initialize on
all the newly created subgrids, this time allocating the data.

I apologize for the incompleteness of this document.. Feel free to
contact David Collins through mailing list if you need to develop a
test using ParallelRootGridIO


