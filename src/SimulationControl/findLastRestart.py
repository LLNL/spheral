import os
import mpi

#-------------------------------------------------------------------------------
# Find the most recent available restart file based on a given Spheral style
# basename.
#-------------------------------------------------------------------------------
def findLastRestart(baseName,
                    rank = mpi.rank,
                    procs = mpi.procs,
                    suffix = ".silo"):

    # Get the list of cycle numbers from the available restart files.
    cycles = findAvailableRestartCycles(baseName, rank, procs, suffix)

    # Finally, we can get the max cycle and return.
    if cycles:
        return max(cycles)
    else:
        return None

#-------------------------------------------------------------------------------
# Find the available set of restart cycles for a given Spheral restart file
# basename.
#-------------------------------------------------------------------------------
def findAvailableRestartCycles(baseName,
                               rank = mpi.rank,
                               procs = mpi.procs,
                               suffix = ".silo"):

    # If we're running mpi we want to make sure we only look at valid restart
    # files for the range we're using.
    if procs > 1:
        baseName += '_rank%i_of_%idomains' % (rank, procs)

    # Get the directory path for the restart files.
    dirPath, fileBase = os.path.split(baseName)
    if dirPath == "":
        dirPath = "."

    # Get a listing of all the files in this directory that match the base name.
    restartFiles = [x for x in os.listdir(dirPath)
                    if x[:len(fileBase)] == fileBase and x.find("_cycle") > 0]

    # Now get the list of cycle numbers from the available restart files.
    cycles = [int(x) for x in
              [(y.split("cycle")[-1]).split(suffix)[0] for y in restartFiles]]
    cycles = list(set(cycles))
    
    # Return that set 'o cycles.
    return cycles

#-------------------------------------------------------------------------------
# Helper method to load up one of the SpheralController auto-generated
# conservation histories from the last restart file under a given path.
#-------------------------------------------------------------------------------
def loadLastConservationHistory(baseName):
    from Spheral import DataBase1d, FlatFileIO
    from SpheralConservation import SpheralConservation

    # Get the last restart cycle under the given base name.
    restartCycle = findLastRestart(baseName)
    if restartCycle is None:
        raise "Unable to find appropriate restart cycle for path %s" % baseName
    print("Selected restart cycle ", restartCycle)

    # Get the directory path for the restart files.
    dirPath, fileBase = os.path.split(baseName)
    if dirPath == "":
        dirPath = "."

    # Find the set of restart files corresponding to the target cycle.
    restartFiles = [x for x in os.listdir(dirPath)
                    if x[:len(fileBase)] == fileBase and x.find("_cycle%i" % restartCycle) > 0]
    assert len(restartFiles) > 0

    # If this is a parallel run, we only need one of these files.
    restartFile = dirPath + "/" + restartFiles[0]

    # Create a FileIO object to read from this file.
    f = FlatFileIO(restartFile, Read)

    # Create a conservation object.
    db = DataBase1d()
    result = SpheralConservation(db, [])

    # Restore the state of the conserve object.
    result.restoreState(f, "restart_cycle%i/control/self/conserve" % restartCycle)
    f.close()

    # That's it.
    return result
