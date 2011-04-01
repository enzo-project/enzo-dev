cimport libc.stdlib

cdef extern from "TopGridData.h":
    struct c_TopGridData "TopGridData":
        Eint   CycleNumber
        Eint   SubcycleNumber
        FLOAT Time
        double CPUTime
        double StartCPUTime
        double LastCycleCPUTime
        Eint ResubmitOn

        # Script names for resubmission to queues and restarting to reduce
        # memory fragmentation.

        char *ResubmitCommand

        # Stopping criteria for TopGrid.

        FLOAT StopTime
        Eint   StopCycle
        Eint   StopSteps
        Eflt StopCPUTime

        # Parameters governing when output is done.

        Eflt TimeLastRestartDump
        Eflt dtRestartDump

        FLOAT TimeLastDataDump
        FLOAT dtDataDump

        FLOAT TimeLastHistoryDump
        FLOAT dtHistoryDump

        FLOAT TimeLastMovieDump
        FLOAT dtMovieDump

        FLOAT TimeLastTracerParticleDump
        FLOAT dtTracerParticleDump

        FLOAT MovieRegionLeftEdge[MAX_DIMENSION]
        FLOAT MovieRegionRightEdge[MAX_DIMENSION]

        FLOAT NewMovieLeftEdge[MAX_DIMENSION]
        FLOAT NewMovieRightEdge[MAX_DIMENSION]

        Eint CycleLastRestartDump
        Eint CycleSkipRestartDump

        Eint CycleLastDataDump
        Eint CycleSkipDataDump

        Eint SubcycleLastDataDump
        Eint SubcycleSkipDataDump

        Eint CycleLastHistoryDump
        Eint CycleSkipHistoryDump
        Eint CycleSkipGlobalDataDump

        Eint OutputFirstTimeAtLevel
        Eint StopFirstTimeAtLevel

        # Parameters governing output names.

        Eint RestartDumpNumber
        Eint DataDumpNumber
        Eint HistoryDumpNumber
        Eint MovieDumpNumber
        Eint TracerParticleDumpNumber

        char *RestartDumpName
        char *DataDumpName
        char *HistoryDumpName
        char *MovieDumpName
        char *TracerParticleDumpName
        char *RedshiftDumpName

        char *RestartDumpDir
        char *DataDumpDir
        char *HistoryDumpDir
        char *MovieDumpDir
        char *TracerParticleDumpDir
        char *RedshiftDumpDir

        char *LocalDir
        char *GlobalDir

        # TopGrid Parameters governing hierarchy

        Eint StaticHierarchy

        # Some grid defining data
        # These are here out of convenience, the real ones are in the grids.

        Eint TopGridRank
        Eint TopGridDims[MAX_DIMENSION]
        #boundary_type  LeftFaceBoundaryCondition[MAX_DIMENSION],
        #              RightFaceBoundaryCondition[MAX_DIMENSION]
        char *BoundaryConditionName

        # Gravity data -- used only for top grid potential field solve

        #gravity_boundary_type GravityBoundary

        # Particle and Particle boundary data. (real one in ExternalBoundary).

        #boundary_type ParticleBoundaryType
        Eint           NumberOfParticles

        Eflt  CourantSafetyNumber
        Eint    PPMFlatteningParameter
        Eint    PPMDiffusionParameter
        Eint    PPMSteepeningParameter

        #AMRHDF5Writer AmiraGrid
        Eint FirstTimestepAfterRestart

cdef class TopGridData:
    cdef c_TopGridData *thisptr
    def __cinit__(self):
        self.thisptr = <c_TopGridData *> libc.stdlib.malloc(sizeof(c_TopGridData))
    def __dealloc__(self):
        libc.stdlib.free(self.thisptr)


    property CycleNumber:
        def __get__(self):
            return self.thisptr.CycleNumber
        def __set__(self, Eint val):
            self.thisptr.CycleNumber = val


    property SubcycleNumber:
        def __get__(self):
            return self.thisptr.SubcycleNumber
        def __set__(self, Eint val):
            self.thisptr.SubcycleNumber = val


    property Time:
        def __get__(self):
            return self.thisptr.Time
        def __set__(self, FLOAT val):
            self.thisptr.Time = val


    property CPUTime:
        def __get__(self):
            return self.thisptr.CPUTime
        def __set__(self, double val):
            self.thisptr.CPUTime = val


    property StartCPUTime:
        def __get__(self):
            return self.thisptr.StartCPUTime
        def __set__(self, double val):
            self.thisptr.StartCPUTime = val


    property LastCycleCPUTime:
        def __get__(self):
            return self.thisptr.LastCycleCPUTime
        def __set__(self, double val):
            self.thisptr.LastCycleCPUTime = val


    property ResubmitOn:
        def __get__(self):
            return self.thisptr.ResubmitOn
        def __set__(self, Eint val):
            self.thisptr.ResubmitOn = val


    property StopTime:
        def __get__(self):
            return self.thisptr.StopTime
        def __set__(self, FLOAT val):
            self.thisptr.StopTime = val


    property StopCycle:
        def __get__(self):
            return self.thisptr.StopCycle
        def __set__(self, Eint val):
            self.thisptr.StopCycle = val


    property StopSteps:
        def __get__(self):
            return self.thisptr.StopSteps
        def __set__(self, Eint val):
            self.thisptr.StopSteps = val


    property StopCPUTime:
        def __get__(self):
            return self.thisptr.StopCPUTime
        def __set__(self, Eflt val):
            self.thisptr.StopCPUTime = val


    property TimeLastRestartDump:
        def __get__(self):
            return self.thisptr.TimeLastRestartDump
        def __set__(self, Eflt val):
            self.thisptr.TimeLastRestartDump = val


    property dtRestartDump:
        def __get__(self):
            return self.thisptr.dtRestartDump
        def __set__(self, Eflt val):
            self.thisptr.dtRestartDump = val


    property TimeLastDataDump:
        def __get__(self):
            return self.thisptr.TimeLastDataDump
        def __set__(self, FLOAT val):
            self.thisptr.TimeLastDataDump = val


    property dtDataDump:
        def __get__(self):
            return self.thisptr.dtDataDump
        def __set__(self, FLOAT val):
            self.thisptr.dtDataDump = val


    property TimeLastHistoryDump:
        def __get__(self):
            return self.thisptr.TimeLastHistoryDump
        def __set__(self, FLOAT val):
            self.thisptr.TimeLastHistoryDump = val


    property dtHistoryDump:
        def __get__(self):
            return self.thisptr.dtHistoryDump
        def __set__(self, FLOAT val):
            self.thisptr.dtHistoryDump = val


    property TimeLastMovieDump:
        def __get__(self):
            return self.thisptr.TimeLastMovieDump
        def __set__(self, FLOAT val):
            self.thisptr.TimeLastMovieDump = val


    property dtMovieDump:
        def __get__(self):
            return self.thisptr.dtMovieDump
        def __set__(self, FLOAT val):
            self.thisptr.dtMovieDump = val


    property TimeLastTracerParticleDump:
        def __get__(self):
            return self.thisptr.TimeLastTracerParticleDump
        def __set__(self, FLOAT val):
            self.thisptr.TimeLastTracerParticleDump = val


    property dtTracerParticleDump:
        def __get__(self):
            return self.thisptr.dtTracerParticleDump
        def __set__(self, FLOAT val):
            self.thisptr.dtTracerParticleDump = val


    property CycleLastRestartDump:
        def __get__(self):
            return self.thisptr.CycleLastRestartDump
        def __set__(self, Eint val):
            self.thisptr.CycleLastRestartDump = val


    property CycleSkipRestartDump:
        def __get__(self):
            return self.thisptr.CycleSkipRestartDump
        def __set__(self, Eint val):
            self.thisptr.CycleSkipRestartDump = val


    property CycleLastDataDump:
        def __get__(self):
            return self.thisptr.CycleLastDataDump
        def __set__(self, Eint val):
            self.thisptr.CycleLastDataDump = val


    property CycleSkipDataDump:
        def __get__(self):
            return self.thisptr.CycleSkipDataDump
        def __set__(self, Eint val):
            self.thisptr.CycleSkipDataDump = val


    property SubcycleLastDataDump:
        def __get__(self):
            return self.thisptr.SubcycleLastDataDump
        def __set__(self, Eint val):
            self.thisptr.SubcycleLastDataDump = val


    property SubcycleSkipDataDump:
        def __get__(self):
            return self.thisptr.SubcycleSkipDataDump
        def __set__(self, Eint val):
            self.thisptr.SubcycleSkipDataDump = val


    property CycleLastHistoryDump:
        def __get__(self):
            return self.thisptr.CycleLastHistoryDump
        def __set__(self, Eint val):
            self.thisptr.CycleLastHistoryDump = val


    property CycleSkipHistoryDump:
        def __get__(self):
            return self.thisptr.CycleSkipHistoryDump
        def __set__(self, Eint val):
            self.thisptr.CycleSkipHistoryDump = val


    property CycleSkipGlobalDataDump:
        def __get__(self):
            return self.thisptr.CycleSkipGlobalDataDump
        def __set__(self, Eint val):
            self.thisptr.CycleSkipGlobalDataDump = val


    property OutputFirstTimeAtLevel:
        def __get__(self):
            return self.thisptr.OutputFirstTimeAtLevel
        def __set__(self, Eint val):
            self.thisptr.OutputFirstTimeAtLevel = val


    property StopFirstTimeAtLevel:
        def __get__(self):
            return self.thisptr.StopFirstTimeAtLevel
        def __set__(self, Eint val):
            self.thisptr.StopFirstTimeAtLevel = val


    property RestartDumpNumber:
        def __get__(self):
            return self.thisptr.RestartDumpNumber
        def __set__(self, Eint val):
            self.thisptr.RestartDumpNumber = val


    property DataDumpNumber:
        def __get__(self):
            return self.thisptr.DataDumpNumber
        def __set__(self, Eint val):
            self.thisptr.DataDumpNumber = val


    property HistoryDumpNumber:
        def __get__(self):
            return self.thisptr.HistoryDumpNumber
        def __set__(self, Eint val):
            self.thisptr.HistoryDumpNumber = val


    property MovieDumpNumber:
        def __get__(self):
            return self.thisptr.MovieDumpNumber
        def __set__(self, Eint val):
            self.thisptr.MovieDumpNumber = val


    property TracerParticleDumpNumber:
        def __get__(self):
            return self.thisptr.TracerParticleDumpNumber
        def __set__(self, Eint val):
            self.thisptr.TracerParticleDumpNumber = val


    property StaticHierarchy:
        def __get__(self):
            return self.thisptr.StaticHierarchy
        def __set__(self, Eint val):
            self.thisptr.StaticHierarchy = val


    property TopGridRank:
        def __get__(self):
            return self.thisptr.TopGridRank
        def __set__(self, Eint val):
            self.thisptr.TopGridRank = val


    property NumberOfParticles:
        def __get__(self):
            return self.thisptr.NumberOfParticles
        def __set__(self, Eint val):
            self.thisptr.NumberOfParticles = val


    property CourantSafetyNumber:
        def __get__(self):
            return self.thisptr.CourantSafetyNumber
        def __set__(self, Eflt val):
            self.thisptr.CourantSafetyNumber = val


    property PPMFlatteningParameter:
        def __get__(self):
            return self.thisptr.PPMFlatteningParameter
        def __set__(self, Eint val):
            self.thisptr.PPMFlatteningParameter = val


    property PPMDiffusionParameter:
        def __get__(self):
            return self.thisptr.PPMDiffusionParameter
        def __set__(self, Eint val):
            self.thisptr.PPMDiffusionParameter = val


    property PPMSteepeningParameter:
        def __get__(self):
            return self.thisptr.PPMSteepeningParameter
        def __set__(self, Eint val):
            self.thisptr.PPMSteepeningParameter = val


    property FirstTimestepAfterRestart:
        def __get__(self):
            return self.thisptr.FirstTimestepAfterRestart
        def __set__(self, Eint val):
            self.thisptr.FirstTimestepAfterRestart = val

    property TopGridDims:
        def __get__(self):
            cdef int i
            a = []
            for i in range(self.thisptr.TopGridRank):
                a.append(self.thisptr.TopGridDims[i])
            return a
        def __set__(self, v):
            cdef int i
            for i in range(self.thisptr.TopGridRank):
                self.thisptr.TopGridDims[i] = v[i]

# SKIPPED: 'FLOAT MovieRegionLeftEdge[MAX_DIMENSION]'
# SKIPPED: 'FLOAT MovieRegionRightEdge[MAX_DIMENSION]'
# SKIPPED: 'FLOAT NewMovieLeftEdge[MAX_DIMENSION]'
# SKIPPED: 'FLOAT NewMovieRightEdge[MAX_DIMENSION]'
# SKIPPED: '#boundary_type  LeftFaceBoundaryCondition[MAX_DIMENSION],'
# SKIPPED: '#              RightFaceBoundaryCondition[MAX_DIMENSION]'
