import mpi
from Spheral import computeVoronoiVolume, HydroFieldNames

def LlyodsAlgorithm(db,hydro,inDomain,rhoAtmo,nPerh,nRelaxIters):
#===============================================================================
# inputs:
#   db ---------- database containing nodeLists for mass correction
#   hydro ------- hydro obj used to access the boundary conditions
#   inDomain ---- logic function returns true if point is inside domain false otherwise
#   rhoAtmo ----- function specifying the density profile
#   nPerh ------- nodes per smoothing length (implementation hack)
#   nRelaxIters - number of centroidal relaxation steps (0 for simple voronoi volume correction)
#===============================================================================
    
    # multiplication factor for centroid relaxation delta
    assert db.nDim in (2,3)
    
    # geometry specific functions
    exec('from Spheral{0}d import SymTensor, IntFieldList, Vector, vector_of_Vector, FacetedVolumeFieldList, vector_of_CellFaceFlagFieldList, vector_of_FacetedVolume, vector_of_Boundary, ScalarFieldList, SymTensorFieldList, vector_of_vector_of_FacetedVolume, vector_of_FacetedVolume'.format(db.nDim))
   
    # Fields we need to create for the voronoi volumme calculation function
    faceted_bounds = vector_of_FacetedVolume()
    bounds = vector_of_Boundary()                      # Bounds enforced by inDomain function
    weight = ScalarFieldList()                         # No weights
    damage = SymTensorFieldList()                      # No damage
    holes = vector_of_vector_of_FacetedVolume([vector_of_FacetedVolume()])
    surfacePoint = IntFieldList() # db.newFluidIntFieldList(0, HydroFieldNames.surfacePoint)
    vol = db.newFluidScalarFieldList(0.0, HydroFieldNames.volume)
    deltaMedian = db.newFluidVectorFieldList(Vector.zero, "centroidal delta")
    etaVoidPoints = db.newFluidvector_of_VectorFieldList(vector_of_Vector(), "eta void points")
    cells = FacetedVolumeFieldList()#db.newFluidFacetedVolumeFieldList(FacetedVolume(), "cells")
    cellFaceFlags = vector_of_CellFaceFlagFieldList()#db.newFluidvector_of_CellFaceFlagFieldList(vector_of_int(), "face flags")
    
    

    # get the bcs from hydro and add in parallel domain bcs if appropriate
    bcs = hydro.boundaryConditions()
    
    # loop through nodeLists and zero out nPerh to 
    # circumvent computeVoronoiVolumes void point feature
    fluidNodeLists = db.fluidNodeLists() 
    #for nodes in fluidNodeLists:

    print(db.numFluidNodeLists,db.numSolidNodeLists)
    for j in range(nRelaxIters): # centroidal relax n times 
        print('centroidal relaxation: iteration {0}'.format(j))

        # centroidal relaxation
        #print('  relax...')
        for k in range(db.numFluidNodeLists):
            for i in range(fluidNodeLists[k].numInternalNodes):
                if inDomain(db.fluidPosition(k,i)): # I used this to handle constant BC
                    db.fluidPosition[k][i] = db.fluidPosition(k,i)+deltaMedian(k,i) # move to centroid
        
        for nodes in fluidNodeLists:
            nodes.numGhostNodes=0
        
        #print('  bcs...')   
        for bc in bcs:
            bc.setAllGhostNodes(db) # generate all the ghost nodes
            bc.finalizeGhostBoundary() # hit the finalize for bcs that need that
        
        #print('  connectivity...')
        db.reinitializeNeighbors() # update neighbors
        db.updateConnectivityMap(False,False) # update connectivity
                #(False,False)  
                    # 1) dont compute connectivity for ghosts 
                    # 2) dont compute overlap connectivity
        
        #print('  volume...')
        computeVoronoiVolume(db.fluidPosition, 
                             db.fluidHfield,
                             db.connectivityMap(), 
                             damage, # believe this is inactive?
                             faceted_bounds, # specify domain extent
                             holes, # polygonal/polyhedral holes in region
                             bounds, # boundaries
                             weight, # weighting function
                             surfacePoint, # between substances
                             vol, # those volumes we want
                             deltaMedian, # distance from centroid
                             etaVoidPoints, # void nodes
                             cells, # polygon/polyhedron object
                             cellFaceFlags)
    for k in range(db.numFluidNodeLists):
        fluidNodeLists[k].nodesPerSmoothingScale = nPerh


    