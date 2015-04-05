#-------------------------------------------------------------------------------
# Dump Spheral meshed data to a set of Silo files obeying the Overlink 
# specification.
#-------------------------------------------------------------------------------
import os, mpi
from Spheral import *
from SpheralModules import silo
from SpheralModules.silo import SiloAttributes as SA

# Parallel info.
domainID = mpi.rank
numDomains = mpi.procs

#-------------------------------------------------------------------------------
# A few integer traits by name that Overlink wants.
#-------------------------------------------------------------------------------
OverlinkAttrs = {"ATTR_NODAL"        : 0,
                 "ATTR_ZONAL"        : 1,
                 "ATTR_FACE"         : 2,
                 "ATTR_EDGE"         : 3,
                 "ATTR_INTENSIVE"    : 0,
                 "ATTR_EXTENSIVE"    : 1,
                 "ATTR_FIRST_ORDER"  : 0,
                 "ATTR_SECOND_ORDER" : 1,
                 "ATTR_INTEGER"      : 0,
                 "ATTR_FLOAT"        : 1,
                 }

#-------------------------------------------------------------------------------
# writeSiloMesh -- this is the one the user should actually call!
#-------------------------------------------------------------------------------
def siloMeshDump(dirName, mesh,
                 index2zone = None,
                 nodeLists = [],
                 label = "Spheral++ generated mesh",
                 time = 0.0,
                 cycle = 0,
                 scalarFields = [],
                 vectorFields = [],
                 tensorFields = [],
                 symTensorFields = [],
                 nodeArrays = None,
                 zoneArrays = None,
                 faceArrays = None,
                 pretendRZ = False):

    assert (isinstance(mesh, polytope.Tessellation2d) or
            isinstance(mesh, polytope.Tessellation3d))
    if isinstance(mesh, polytope.Tessellation2d):
        nDim = 2
    else:
        nDim = 3

    # We can only pretend this is an RZ mesh if it's 2D.
    assert (not pretendRZ) or (nDim == 2)
    
    # Make sure the requested directory exists!
    if domainID == 0:
        if not os.path.exists(dirName):
            os.makedirs(dirName)
    mpi.barrier()

    # Extract all the fields we're going to write.
    fieldwad = extractFieldComponents(nodeLists, time, cycle,
                                      scalarFields, vectorFields, tensorFields, symTensorFields)

    # If we have index2zone, remove any redundant values.
    if index2zone:
        ntot = sum([nodes.numInternalNodes for nodes in nodeLists])
        assert len(index2zone) == ntot
        ntarget = len(mesh.cells)
        assert max(index2zone) + 1 == ntarget
        for name, desc, ftype, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
            newsubvars = []
            for subname, vals in subvars:
                valcopy = list(vals)
                vals.resize(ntarget)
                for i, vali in enumerate(valcopy):
                    vals[index2zone[i]] = vali
                assert len(vals) == ntarget
        #         newsubvars.append((subname, newvals))
        #     assert len(newsubvars) == len(subvars)
        #     newfieldwad.append([name, desc, ftype, optlistDef, optlistMV, optlistVar, newsubvars])
        # assert len(newfieldwad) == len(fieldwad)
        # fieldwad = newfieldwad
    else:
        index2zone = range(len(mesh.cells))

    # If we're domain 0 we write the master file.
    masterfile = None
    if domainID == 0:
        masterfile = writeMasterMeshSiloFile(dirName, label, nodeLists, time, cycle, fieldwad)
    masterfile = mpi.bcast(masterfile)

    # Each domain writes it's domain file.
    writeDomainMeshSiloFile(dirName, mesh, index2zone, label, nodeLists, time, cycle, fieldwad,
                            pretendRZ, nodeArrays, zoneArrays, faceArrays)


    # That's it.
    return masterfile

#-------------------------------------------------------------------------------
# Extract the fields we're going to write.  This requires exploding vector and
# tensor fields into their components.
#-------------------------------------------------------------------------------
def extractFieldComponents(nodeLists, time, cycle,
                           scalarFields, vectorFields, tensorFields, symTensorFields):
    result = extractFields(nodeLists, time, cycle, scalarFields,
                           metaDataScalarField, extractScalarField, dummyScalarField)
    result += extractFields(nodeLists, time, cycle, vectorFields,
                            metaDataVectorField, extractVectorField, dummyVectorField)
    result += extractFields(nodeLists, time, cycle, tensorFields,
                            metaDataTensorField, extractTensorField, dummyTensorField)
    result += extractFields(nodeLists, time, cycle, symTensorFields,
                            metaDataTensorField, extractTensorField, dummyTensorField)
    return result

#-------------------------------------------------------------------------------
# Write the master file.
#-------------------------------------------------------------------------------
def writeMasterMeshSiloFile(dirName, label, nodeLists, time, cycle, fieldwad,
                            meshType = SA._DB_UCDMESH):

    nullOpts = silo.DBoptlist()

    # Create the master file.
    p0, p1 = os.path.split(dirName)
    fileName = os.path.join(dirName, "OvlTop.silo")
    db = silo.DBCreate(fileName, 
                       SA._DB_CLOBBER, SA._DB_LOCAL, label, SA._DB_HDF5)

    # Pattern for constructing per domain variables.
    domainNamePatterns = [("%s/domain%i.silo:" % (p1, i)) + "%s" for i in xrange(numDomains)]

    # Write the domain file names and types.
    domainNames = vector_of_string()
    meshTypes = vector_of_int()
    i = 0
    for p in domainNamePatterns:
        domainNames.append(p % "MESH")
        meshTypes.append(meshType)
    optlist = silo.DBoptlist(1024)
    assert optlist.addOption(SA._DBOPT_CYCLE, cycle) == 0
    assert optlist.addOption(SA._DBOPT_DTIME, time) == 0
    assert silo.DBPutMultimesh(db, "MMESH", domainNames, meshTypes, optlist) == 0

    # Extract the material names, and write per material info if any.
    if len(nodeLists) > 0:

        # Write material names.
        materialNames = [x.name for x in nodeLists]
        material_names = vector_of_string()
        matnames = vector_of_string()
        matnos = vector_of_int()
        for p in domainNamePatterns:
            material_names.append(p % "MATERIAL")
        for (name, i) in zip(materialNames, range(len(materialNames))):
            matnames.append(name)
            matnos.append(i)
        assert len(material_names) == numDomains
        assert len(matnames) == len(nodeLists)
        assert len(matnos) == len(nodeLists)
        optlist = silo.DBoptlist(1024)
        assert optlist.addOption(SA._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(SA._DBOPT_DTIME, time) == 0
        assert optlist.addOption(SA._DBOPT_MATNAMES, SA._DBOPT_NMATNOS, matnames) == 0
        assert optlist.addOption(SA._DBOPT_MATNOS, SA._DBOPT_NMATNOS, matnos) == 0
        assert silo.DBPutMultimat(db, "MMATERIAL", material_names, optlist) == 0

        # Write the variable descriptions for non-scalar variables (vector and tensors).
        writeDefvars(db, fieldwad)

        # Write the attributes for each variable.
        # This is largely an Overlink thing, which only needs the components of types like vectors.
        if len(fieldwad) > 0:
            attrNames = vector_of_string()
            attrs = vector_of_vector_of_int()
            thpt = vector_of_int()
            thpt.append(OverlinkAttrs["ATTR_ZONAL"])
            thpt.append(OverlinkAttrs["ATTR_INTENSIVE"])
            thpt.append(OverlinkAttrs["ATTR_SECOND_ORDER"])
            thpt.append(0)
            thpt.append(OverlinkAttrs["ATTR_FLOAT"])
            for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
                for subname, vals in subvars:
                    attrNames.append(subname)
                    attrs.append(thpt)
            assert len(attrNames) == len(attrs)
            silo.DBPutCompoundarray(db, "VAR_ATTRIBUTES", attrNames, attrs, nullOpts)

        # Write the variables descriptors.
        ucdTypes = vector_of_int(numDomains, SA._DB_UCDVAR)
        for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
            domainVarNames = vector_of_string()
            for p in domainNamePatterns:
                domainVarNames.append(p % name)
            assert len(domainVarNames) == numDomains
            assert silo.DBPutMultivar(db, name, domainVarNames, ucdTypes, optlistMV) == 0
            if desc != None:
                for subname, vals in subvars:
                    domainVarNames = vector_of_string()
                    for p in domainNamePatterns:
                        domainVarNames.append(p % subname)
                    assert len(domainVarNames) == numDomains
                    assert silo.DBPutMultivar(db, subname, domainVarNames, ucdTypes, optlistVar) == 0

    # Make a convenient symlink for the master file.
    linkfile = p1 + ".silo"
    fullpath = os.path.join(p0, linkfile)
    if os.path.exists(fullpath):
        os.remove(fullpath)
    targetfile = os.path.join(p1, "OvlTop.silo")
    os.symlink(targetfile, os.path.join(p0, linkfile))

    # That's it.
    assert silo.DBClose(db) == 0
    del db
    return linkfile

#-------------------------------------------------------------------------------
# Write the domain file.
#-------------------------------------------------------------------------------
def writeDomainMeshSiloFile(dirName, mesh, index2zone, label, nodeLists, time, cycle, fieldwad,
                            pretendRZ, nodeArrays, zoneArrays, faceArrays,
                            meshType = SA._DB_UCDMESH):

    # Create the file.
    fileName = os.path.join(dirName, "domain%i.silo" % domainID)
    db = silo.DBCreate(fileName, 
                       SA._DB_CLOBBER, SA._DB_LOCAL, label, SA._DB_HDF5)
    nullOpts = silo.DBoptlist()

    # Determine our dimensionality
    if isinstance(mesh, polytope.Tessellation2d):
        nDim = 2
    else:
        assert isinstance(mesh, polytope.Tessellation3d)
        nDim = 3

    # Write a Polygonal zone list.
    zonelistName = { 2 : "zonelist",
                     3 : "PHzonelist" }
    numZones = len(mesh.cells)
    if nDim == 2:

        # Read out the zone nodes.  We rely on these already being arranged
        # counter-clockwise.
        zoneNodes = vector_of_vector_of_int()
        shapesize = vector_of_int()
        for zoneID in xrange(numZones):
            zone = mesh.cells[zoneID]
            nodes = vector_of_int()
            for iface in zone:
                if iface < 0:
                    nodes.append(mesh.faces[~iface][1])
                else:
                    nodes.append(mesh.faces[iface][0])
            zoneNodes.append(nodes)
            shapesize.append(len(nodes))
        assert len(zoneNodes) == numZones
        assert len(shapesize) == numZones

        assert silo.DBPutZonelist2(db, zonelistName[nDim], nDim, zoneNodes, 0, 0,
                                   vector_of_int(numZones, SA._DB_ZONETYPE_POLYGON),
                                   shapesize,
                                   vector_of_int(numZones, 1),
                                   nullOpts) == 0

    # Write a Polyhedral zone list.
    if nDim == 3:

        # Construct the face-node lists.
        numFaces = len(mesh.faces)
        faceNodes = vector_of_vector_of_int(numFaces)
        for iface in xrange(numFaces):
            for j in xrange(len(mesh.faces[iface])):
                faceNodes[iface].append(mesh.faces[iface][j])
            assert len(faceNodes[-1]) == face.numNodes
        assert len(faceNodes) == numFaces

        # Construct the zone-face list.  We use the ones complement of a face ID
        # to indicate that face needs to be reversed in reference to this zone.
        # This is the same convention as polytope, so just copy it.
        zoneFaces = mesh.cells
        assert len(zoneFaces) == numZones

        assert silo.DBPutPHZonelist(db, zonelistName[nDim], faceNodes, zoneFaces, 0, (numZones - 1), nullOpts) == 0

    # Construct the mesh node coordinates.
    assert len(mesh.nodes) % nDim == 0
    numNodes = len(mesh.nodes)/nDim
    coords = vector_of_vector_of_double(nDim, vector_of_double(numNodes))
    for nodeID in xrange(numNodes):
        for idim in xrange(nDim):
            coords[idim][nodeID] = mesh.nodes[nDim*nodeID + idim]
    assert len(coords) == nDim

    # Write the mesh itself.
    meshOpts = silo.DBoptlist(1024)
    assert meshOpts.addOption(SA._DBOPT_CYCLE, cycle) == 0
    assert meshOpts.addOption(SA._DBOPT_DTIME, time) == 0
    assert meshOpts.addOption(SA._DBOPT_COORDSYS, SA._DB_CARTESIAN) == 0
    assert meshOpts.addOption(SA._DBOPT_NSPACE, nDim) == 0
    assert meshOpts.addOption(SA._DBOPT_TV_CONNECTIVITY, 1) == 0
    if nDim == 2:
        if pretendRZ:
            assert meshOpts.addOption(SA._DBOPT_COORDSYS, SA._DB_CYLINDRICAL) == 0
            assert meshOpts.addOption(SA._DBOPT_XLABEL, "z") == 0
            assert meshOpts.addOption(SA._DBOPT_YLABEL, "r") == 0
        assert silo.DBPutUcdmesh(db, "MESH", coords, numZones, zonelistName[nDim], "NULL", meshOpts) == 0
    else:
        assert meshOpts.addOption(SA._DBOPT_PHZONELIST, zonelistName[nDim]) == 0
        assert silo.DBPutUcdmesh(db, "MESH", coords, numZones, "NULL", "NULL", meshOpts) == 0

    # Write materials.
    if len(nodeLists) > 0:
        matnos = vector_of_int()
        for i in xrange(len(nodeLists)):
            matnos.append(i)
        assert len(matnos) == len(nodeLists)
        matlist = vector_of_int(numZones)
        matnames = vector_of_string()
        offset = 0
        for (nodeList, imat) in zip(nodeLists, xrange(len(nodeLists))):
            for i in xrange(nodeList.numInternalNodes):
                matlist[index2zone[offset + i]] = imat
            matnames.append(nodeList.name)
            offset += nodeList.numInternalNodes
        assert len(matlist) == numZones
        assert len(matnames) == len(nodeLists)
        matOpts = silo.DBoptlist(1024)
        assert matOpts.addOption(SA._DBOPT_CYCLE, cycle) == 0
        assert matOpts.addOption(SA._DBOPT_DTIME, time) == 0
        assert matOpts.addOption(SA._DBOPT_MATNAMES, SA._DBOPT_NMATNOS, matnames) == 0
        assert silo.DBPutMaterial(db, "MATERIAL", "MESH", matnos, matlist,
                                  vector_of_int(), vector_of_int(), vector_of_int(), vector_of_double(),
                                  matOpts) == 0

        # Write the variable descriptions for non-scalar variables (vector and tensors).
        writeDefvars(db, fieldwad)

        # Write the field components.
        centering = SA._DB_ZONECENT
        varOpts = silo.DBoptlist(1024)
        assert varOpts.addOption(SA._DBOPT_CYCLE, cycle) == 0
        assert varOpts.addOption(SA._DBOPT_DTIME, time) == 0
        for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
            for subname, vals in subvars:
                if len(vals) > 0:
                    assert silo.DBPutUcdvar1(db, subname, "MESH", vals, vector_of_double(), centering, varOpts) == 0

    # Write the set of neighbor domains.
    thpt = vector_of_vector_of_int()
    thpt.append(vector_of_int(1, len(mesh.neighborDomains)))
    thpt.append(vector_of_int())
    for i in mesh.neighborDomains:
        thpt[-1].append(i)
    elemNames = vector_of_string()
    elemNames.append("num neighbor domains")
    elemNames.append("neighbor domains")
    assert silo.DBPutCompoundarray(db, "DOMAIN_NEIGHBOR_NUMS", elemNames, thpt, nullOpts) == 0

    # Write the shared nodes for each neighbor domain.
    sharedNodes = mesh.sharedNodes
    for ineighborDomain in xrange(len(mesh.neighborDomains)):
        nodes = vector_of_int()
        for i in sharedNodes[ineighborDomain]:
            nodes.append(i)
        assert len(nodes) == len(sharedNodes[ineighborDomain])
        assert silo.DBPutCompoundarray(db, "DOMAIN_NEIGHBOR%i" % ineighborDomain,
                                       vector_of_string(1, "shared_nodes"),
                                       vector_of_vector_of_int(1, nodes),
                                       nullOpts) == 0

    # If requested, write out annotations for the nodes, zones, and faces.
    if (not (nodeArrays is None) or
        not (zoneArrays is None) or
        not (faceArrays is None)):
        names = vector_of_string()
        values = vector_of_vector_of_int()
        if not (nodeArrays is None):
            for pair in nodeArrays:
                assert len(pair) == 2
                if len(pair[1]) > 0:
                    names.append(pair[0] + "_node")
                    values.append(vector_of_int())
                    for i in pair[1]:
                        values[-1].append(i)
                    assert len(values[-1]) == len(pair[1])
        if not (zoneArrays is None):
            for pair in zoneArrays:
                assert len(pair) == 2
                if len(pair[1]) > 0:
                    names.append(pair[0] + "_zone")
                    values.append(vector_of_int())
                    for i in pair[1]:
                        values[-1].append(i)
                    assert len(values[-1]) == len(pair[1])
        if not (faceArrays is None):
            for pair in faceArrays:
                assert len(pair) == 2
                if len(pair[1]) > 0:
                    names.append(pair[0] + "_face")
                    values.append(vector_of_int())
                    for i in pair[1]:
                        values[-1].append(i)
                    assert len(values[-1]) == len(pair[1])
        assert len(names) == len(values)
        if len(names) > 0:
            assert silo.DBPutCompoundarray(db, "ANNOTATION_INT", names, values, nullOpts) == 0

    # That's it.
    assert silo.DBClose(db) == 0
    del db
    return

#-------------------------------------------------------------------------------
# Generic master method to extract full sets of field values.
#-------------------------------------------------------------------------------
def extractFields(nodeLists, time, cycle, fields,
                  metaDataMethod,
                  extractMethod,
                  dudMethod):
    result = []
    if len(fields) > 0:
        dim = dimension(fields[0])

        # Figure out how many values we should have for each field.
        nvals = sum([nodes.numInternalNodes for nodes in nodeLists])

        # Group the fields by name (across all NodeLists).
        fieldsByName = {}
        for f in fields:
            if f.name in fieldsByName:
                fieldsByName[f.name][f.nodeList().name] = f
            else:
                fieldsByName[f.name] = {f.nodeList().name : f}

        # Go through and build up the full set of values.
        for fieldName in fieldsByName:
            subfields = fieldsByName[fieldName]
            assert len(subfields) <= len(nodeLists)

            # Build the meta-data for this field entry.
            name = str(fieldName).replace(" ", "_")
            varDef, varType, optlistDef, optlistMV, optlistVar = metaDataMethod(name, time, cycle, dim)

            # Build the complete values for this field.
            vals = []
            for nodes in nodeLists:
                if nodes.name in subfields:
                    vals = extractMethod(name, subfields[nodes.name].internalValues(), vals, dim)
                else:
                    vals = dudMethod(name, nodes.numInternalNodes, vals, dim)
            for thpt in vals:
                assert len(thpt) == 2
                assert len(thpt[1]) == nvals

            result.append((name, varDef, varType, optlistDef, optlistMV, optlistVar, vals))

    return result

#-------------------------------------------------------------------------------
# Scalar field components.
#-------------------------------------------------------------------------------
def extractScalarField(name, field, vals, dim):
    if vals == []:
        vals = [[name, field]]
    else:
        vals[0][1] += field
    return vals

def dummyScalarField(name, n, vals, dim):
    if vals == []:
        vals = [[name, vector_of_double(n, 0.0)]]
    else:
        vals[0][1] += vector_of_double(n, 0.0)
    return vals

def metaDataScalarField(name, time, cycle, dim):
    optlistDef = None
    optlistMV = silo.DBoptlist()
    optlistVar = silo.DBoptlist()
    for optlist in (optlistMV, optlistVar):
        assert optlist.addOption(SA._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(SA._DBOPT_DTIME, time) == 0
        assert optlist.addOption(SA._DBOPT_TENSOR_RANK, SA._DB_VARTYPE_SCALAR) == 0
    return (None, SA._DB_VARTYPE_SCALAR, optlistDef, optlistMV, optlistVar)

#-------------------------------------------------------------------------------
# Vector field components.
#-------------------------------------------------------------------------------
def extractVectorField(name, field, vals, dim):
    assert len(vals) == dim or len(vals) == 0
    assert dim in (2,3)

    if dim == 2:
        if vals == []:
            vals = [["%s_x" % name, vector_of_double()],
                    ["%s_y" % name, vector_of_double()]]
        for v in field:
            vals[0][1].append(v.x)
            vals[1][1].append(v.y)

    else:
        if vals == []:
            vals = [["%s_x" % name, vector_of_double()],
                    ["%s_y" % name, vector_of_double()],
                    ["%s_z" % name, vector_of_double()]]
        for v in field:
            vals[0][1].append(v.x)
            vals[1][1].append(v.y)
            vals[2][1].append(v.z)

    return vals

def dummyVectorField(name, n, vals, dim):
    assert len(vals) == dim or len(vals) == 0
    assert dim in (2,3)
    if vals == []:
        if dim == 2:
            vals = [["%s_x" % name, vector_of_double(n, 0.0)],
                    ["%s_y" % name, vector_of_double(n, 0.0)]]
        else:
            vals = [["%s_x" % name, vector_of_double(n, 0.0)],
                    ["%s_y" % name, vector_of_double(n, 0.0)],
                    ["%s_z" % name, vector_of_double(n, 0.0)]]
    else:
        for i in xrange(dim):
            vals[i][1] += vector_of_double(n, 0.0)
    return vals

def metaDataVectorField(name, time, cycle, dim):
    assert dim in (2,3)
    optlistDef = silo.DBoptlist()
    optlistMV = silo.DBoptlist()
    optlistVar = silo.DBoptlist()
    for optlist in (optlistDef, optlistMV, optlistVar):
        assert optlist.addOption(SA._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(SA._DBOPT_DTIME, time) == 0
    assert optlistMV.addOption(SA._DBOPT_TENSOR_RANK, SA._DB_VARTYPE_VECTOR) == 0
    assert optlistVar.addOption(SA._DBOPT_HIDE_FROM_GUI, 1) == 0
    assert optlistVar.addOption(SA._DBOPT_TENSOR_RANK, SA._DB_VARTYPE_SCALAR) == 0

    if dim == 2:
        return ("{%s_x, %s_y}" % (name, name), SA._DB_VARTYPE_VECTOR,
                optlistDef, optlistMV, optlistVar)
    else:
        return ("{%s_x, %s_y, %s_z}" % (name, name, name), SA._DB_VARTYPE_VECTOR,
                optlistDef, optlistMV, optlistVar)

#-------------------------------------------------------------------------------
# Tensor fields components.
#-------------------------------------------------------------------------------
def extractTensorField(name, field, vals, dim):
    assert len(vals) == dim*dim or len(vals) == 0
    assert dim in (2,3)
    if dim == 2:
        if vals == []:
            vals = [["%s_xx" % name, vector_of_double()],
                    ["%s_xy" % name, vector_of_double()],
                    ["%s_yx" % name, vector_of_double()],
                    ["%s_yy" % name, vector_of_double()]]
        for t in field:
            vals[0][1].append(t.xx); vals[1][1].append(t.xy)
            vals[2][1].append(t.yx); vals[3][1].append(t.yy)

    else:
        if vals == []:
            vals = [["%s_xx" % name, vector_of_double()],
                    ["%s_xy" % name, vector_of_double()],
                    ["%s_xz" % name, vector_of_double()],
                    ["%s_yx" % name, vector_of_double()],
                    ["%s_yy" % name, vector_of_double()],
                    ["%s_yz" % name, vector_of_double()],
                    ["%s_zx" % name, vector_of_double()],
                    ["%s_zy" % name, vector_of_double()],
                    ["%s_zz" % name, vector_of_double()]]
        for t in field:
            vals[0][1].append(t.xx); vals[1][1].append(t.xy); vals[2][1].append(t.xz)
            vals[3][1].append(t.yx); vals[4][1].append(t.yy); vals[5][1].append(t.yz)
            vals[6][1].append(t.zx); vals[7][1].append(t.zy); vals[8][1].append(t.zz)

    return vals

def dummyTensorField(name, n, vals, dim):
    assert len(vals) == dim*dim or len(vals) == 0
    assert dim in (2,3)
    if vals == []:
        if dim == 2:
            vals = [["%s_xx" % name, vector_of_double(n, 0.0)],
                    ["%s_xy" % name, vector_of_double(n, 0.0)],
                    ["%s_yx" % name, vector_of_double(n, 0.0)],
                    ["%s_yy" % name, vector_of_double(n, 0.0)]]
        else:
            vals = [["%s_xx" % name, vector_of_double(n, 0.0)],
                    ["%s_xy" % name, vector_of_double(n, 0.0)],
                    ["%s_xz" % name, vector_of_double(n, 0.0)],
                    ["%s_yx" % name, vector_of_double(n, 0.0)],
                    ["%s_yy" % name, vector_of_double(n, 0.0)],
                    ["%s_yz" % name, vector_of_double(n, 0.0)],
                    ["%s_zx" % name, vector_of_double(n, 0.0)],
                    ["%s_zy" % name, vector_of_double(n, 0.0)],
                    ["%s_zz" % name, vector_of_double(n, 0.0)]]
    else:
        for i in xrange(dim*dim):
            vals[i][1] += vector_of_double(n, 0.0)
    return vals

def metaDataTensorField(name, time, cycle, dim):
    assert dim in (2,3)
    optlistDef = silo.DBoptlist()
    optlistMV = silo.DBoptlist()
    optlistVar = silo.DBoptlist()
    for optlist in (optlistDef, optlistMV, optlistVar):
        assert optlist.addOption(SA._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(SA._DBOPT_DTIME, time) == 0
    assert optlistMV.addOption(SA._DBOPT_TENSOR_RANK, SA._DB_VARTYPE_TENSOR) == 0
    assert optlistVar.addOption(SA._DBOPT_HIDE_FROM_GUI, 1) == 0
    assert optlistVar.addOption(SA._DBOPT_TENSOR_RANK, SA._DB_VARTYPE_SCALAR) == 0

    if dim == 2:
        return ("{{%s_xx, %s_xy}, {%s_yx, %s_yy}}" % (name, name, name, name), SA._DB_VARTYPE_TENSOR,
                optlistDef, optlistMV, optlistVar)
    else:
        return ("{{%s_xx, %s_xy, %s_xz}, {%s_yx, %s_yy, %s_yz}, {%s_zx, %s_zy, %s_zz}}" % (name, name, name,
                                                                                           name, name, name,
                                                                                           name, name, name),
                SA._DB_VARTYPE_TENSOR, optlistDef, optlistMV, optlistVar)

#-------------------------------------------------------------------------------
# Write the variable descriptors for non-scalar types (vector and tensor).
#-------------------------------------------------------------------------------
def writeDefvars(db, fieldwad):
    names, defs, types, opts = vector_of_string(), vector_of_string(), vector_of_int(), vector_of_DBoptlist()
    for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
        if desc != None:
            assert optlistDef != None
            assert len(subvars) > 1
            names.append(name)
            defs.append(desc)
            types.append(type)
            opts.append(optlistDef)
    if len(names) > 0:
        assert silo.DBPutDefvars(db, "VARDEFS", names, types, defs, opts) == 0

#-------------------------------------------------------------------------------
# Extract the dimensionality of a field.
#-------------------------------------------------------------------------------
def dimension(field):
    if "Field2d" in str(field):
        return 2
    elif "Field3d" in str(field):
        return 3
    else:
        assert False
