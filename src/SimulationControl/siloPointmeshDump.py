#-------------------------------------------------------------------------------
# Dump Spheral point data to a set of Silo files using the Pointmesh silo 
# structures.
#-------------------------------------------------------------------------------
import os, mpi
from Spheral import *
from SpheralModules import silo

#-------------------------------------------------------------------------------
# siloPointMeshDump -- this is the one the user should actually call!
#-------------------------------------------------------------------------------
def siloPointmeshDump(baseName, 
                      fields = [],
                      fieldLists = [],
                      procDirBaseName = "proc-%06i",
                      label = "Spheral++ point mesh",
                      time = 0.0,
                      cycle = 0):

    # You have to give us something!
    if len(fields) + len(fieldLists) == 0:
        raise ValueError, "siloPointmeshDump called with no information to write."

    # Get the set of NodeLists we're working with.
    nodeListsDict = {}
    for field in fields:
        nodes = field.nodeList()
        nodeListsDict[nodes.name] = nodes
    for fl in fieldLists:
        nls = fl.nodeListPtrs()
        for nodes in nls:
            nodeListsDict[nodes.name] = nodes
    nodeLists = [nodeListsDict[name] for name in nodeListsDict]
    assert len(nodeLists) > 0

    # If needed, create the subdirectories for the processor files.
    for i in xrange(mpi.procs):
        dire = procDirBaseName % i
        if not os.path.exists(dire):
            os.makedirs(dire)
    mpi.barrier()

    # We can only pretend this is an RZ mesh if it's 2D.
    if isinstance(nodeLists[0], NodeList2d):
        ndim = 2
    elif isinstance(nodeLists[0], NodeList3d):
        ndim = 3
    else:
        raise ValueError, "You need to provide 2D or 3D information for siloPointMeshDump."
    
    # Characterize the fields we're going to write.
    allfields = fields[:]
    for fl in fieldLists:
        for f in fl:
            allfields.append(f)
    intFields, scalarFields, vectorFields, tensorFields, symTensorFields = [], [], [], [], []
    for f in allfields:
        if isinstance(f, eval("IntField%id" % ndim)):
            intFields.append(f)
        elif isinstance(f, eval("ScalarField%id" % ndim)):
            scalarFields.append(f)
        elif isinstance(f, eval("VectorField%id" % ndim)):
            vectorFields.append(f)
        elif isinstance(f, eval("VectorField%id" % ndim)):
            vectorFields.append(f)
        elif isinstance(f, eval("TensorField%id" % ndim)):
            tensorFields.append(f)
        elif isinstance(f, eval("SymTensorField%id" % ndim)):
            symTensorFields.append(f)
        else:
            print "siloPointmeshDump WARNING: ignoring unknown field type."

    # Extract all the fields we're going to write.
    fieldwad = extractFieldComponents(nodeLists, time, cycle, 
                                      intFields, scalarFields, vectorFields, tensorFields, symTensorFields)

    # If we're domain 0 we write the master file.
    if mpi.rank == 0:
        writeMasterSiloFile(ndim, baseName, procDirBaseName, nodeLists, label, time, cycle, fieldwad)

    # Each domain writes it's domain file.
    writeDomainSiloFile(ndim, baseName, procDirBaseName, nodeLists, label, time, cycle, fieldwad)

#-------------------------------------------------------------------------------
# Extract the fields we're going to write.  This requires exploding vector and
# tensor fields into their components.
#-------------------------------------------------------------------------------
def extractFieldComponents(nodeLists, time, cycle,
                           intFields, scalarFields, vectorFields, tensorFields, symTensorFields):
    result = extractFields(nodeLists, time, cycle, intFields,
                           metaDataIntField, extractIntField, dummyIntField)
    result += extractFields(nodeLists, time, cycle, scalarFields,
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
def writeMasterSiloFile(ndim, baseName, procDirBaseName, nodeLists,
                        label, time, cycle, fieldwad):

    nullOpts = silo.DBoptlist()

    # Create the master file.
    fileName = baseName + ".silo"
    db = silo.DBCreate(fileName, 
                       silo._DB_CLOBBER, silo._DB_LOCAL, label, silo._DB_HDF5)

    # Pattern for constructing per domain variables.
    domainNamePatterns = [os.path.join(procDirBaseName % i, baseName + ".silo:%s") for i in xrange(mpi.procs)]

    # Write the domain file names and types.
    domainNames = vector_of_string()
    meshTypes = vector_of_int(mpi.procs, silo._DB_POINTMESH)
    for p in domainNamePatterns:
        domainNames.append(p % "MESH")
    optlist = silo.DBoptlist(1024)
    assert optlist.addOption(silo._DBOPT_CYCLE, cycle) == 0
    assert optlist.addOption(silo._DBOPT_DTIME, time) == 0
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
        assert len(material_names) == mpi.procs
        assert len(matnames) == len(nodeLists)
        assert len(matnos) == len(nodeLists)
        optlist = silo.DBoptlist(1024)
        assert optlist.addOption(silo._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo._DBOPT_DTIME, time) == 0
        assert optlist.addOption(silo._DBOPT_MATNAMES, silo._DBOPT_NMATNOS, matnames) == 0
        assert optlist.addOption(silo._DBOPT_MATNOS, silo._DBOPT_NMATNOS, matnos) == 0
        assert silo.DBPutMultimat(db, "MMATERIAL", material_names, optlist) == 0

        # Write the variable descriptions for non-scalar variables (vector and tensors).
        writeDefvars(db, fieldwad)

        # Write the variables descriptors.
        ucdTypes = vector_of_int(mpi.procs, silo._DB_UCDVAR)
        for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
            domainVarNames = vector_of_string()
            for p in domainNamePatterns:
                domainVarNames.append(p % name)
            assert len(domainVarNames) == mpi.procs
            assert silo.DBPutMultivar(db, name, domainVarNames, ucdTypes, optlistMV) == 0
            if desc != None:
                for subname, vals in subvars:
                    domainVarNames = vector_of_string()
                    for p in domainNamePatterns:
                        domainVarNames.append(p % subname)
                    assert len(domainVarNames) == mpi.procs
                    assert silo.DBPutMultivar(db, subname, domainVarNames, ucdTypes, optlistVar) == 0

    # That's it.
    assert silo.DBClose(db) == 0
    del db
    return

#-------------------------------------------------------------------------------
# Write the domain file.
#-------------------------------------------------------------------------------
def writeDomainSiloFile(ndim, baseName, procDirBaseName, nodeLists, label, time, cycle, fieldwad):

    # Create the file.
    fileName = os.path.join(procDirBaseName % mpi.rank, baseName + ".silo")
    db = silo.DBCreate(fileName, 
                       silo._DB_CLOBBER, silo._DB_LOCAL, label, silo._DB_HDF5)
    nullOpts = silo.DBoptlist()
    print "db is ", db

    # Read the per material info.
    ntot = sum([n.numInternalNodes for n in nodeLists])
    coords = vector_of_vector_of_double(ndim)
    for nodes in nodeLists:
        pos = nodes.positions().internalValues()
        n = len(pos)
        for j in xrange(ndim):
            for i in xrange(n):
                coords[j].append(pos[i][j])
    for j in xrange(ndim):
        assert len(coords[j]) == ntot

    # Write the Pointmesh.
    meshOpts = silo.DBoptlist(1024)
    assert meshOpts.addOption(silo._DBOPT_CYCLE, cycle) == 0
    assert meshOpts.addOption(silo._DBOPT_DTIME, time) == 0
    assert silo.DBPutPointmesh(db, "MESH", coords, meshOpts) == 0

    # Write materials.
    matnos = vector_of_int()
    for i in xrange(len(nodeLists)):
        matnos.append(i)
    assert len(matnos) == len(nodeLists)
    matlist = vector_of_int()
    matnames = vector_of_string()
    for (nodeList, imat) in zip(nodeLists, xrange(len(nodeLists))):
        matlist += vector_of_int(nodeList.numInternalNodes, imat)
        matnames.append(nodeList.name)
    assert len(matlist) == ntot
    assert len(matnames) == len(nodeLists)
    matOpts = silo.DBoptlist(1024)
    assert matOpts.addOption(silo._DBOPT_CYCLE, cycle) == 0
    assert matOpts.addOption(silo._DBOPT_DTIME, time) == 0
    assert matOpts.addOption(silo._DBOPT_MATNAMES, silo._DBOPT_NMATNOS, matnames) == 0
    assert silo.DBPutMaterial(db, "MATERIAL", "MESH", matnos, matlist,
                              vector_of_int(), vector_of_int(), vector_of_int(), vector_of_double(),
                              matOpts) == 0

    # Write the variable descriptions for non-scalar variables (vector and tensors).
    writeDefvars(db, fieldwad)

    # Write the field components.
    centering = silo._DB_ZONECENT
    varOpts = silo.DBoptlist(1024)
    assert varOpts.addOption(silo._DBOPT_CYCLE, cycle) == 0
    assert varOpts.addOption(silo._DBOPT_DTIME, time) == 0
    for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
        for subname, vals in subvars:
            if len(vals) > 0:
                assert silo.DBPutPointvar1(db, subname, "MESH", vals, varOpts) == 0

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
# Int field components.
#-------------------------------------------------------------------------------
def extractIntField(name, field, vals, dim):
    if vals == []:
        vals = [[name, field]]
    else:
        vals[0][1] += field
    return vals

def dummyIntField(name, n, vals, dim):
    if vals == []:
        vals = [[name, vector_of_int(n, 0)]]
    else:
        vals[0][1] += vector_of_int(n, 0)
    return vals

def metaDataIntField(name, time, cycle, dim):
    optlistDef = None
    optlistMV = silo.DBoptlist()
    optlistVar = silo.DBoptlist()
    for optlist in (optlistMV, optlistVar):
        assert optlist.addOption(silo._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo._DBOPT_DTIME, time) == 0
        assert optlist.addOption(silo._DBOPT_TENSOR_RANK, silo._DB_VARTYPE_SCALAR) == 0
    return (None, silo._DB_VARTYPE_INT, optlistDef, optlistMV, optlistVar)

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
        assert optlist.addOption(silo._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo._DBOPT_DTIME, time) == 0
        assert optlist.addOption(silo._DBOPT_TENSOR_RANK, silo._DB_VARTYPE_SCALAR) == 0
    return (None, silo._DB_VARTYPE_SCALAR, optlistDef, optlistMV, optlistVar)

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
        assert optlist.addOption(silo._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo._DBOPT_DTIME, time) == 0
    assert optlistMV.addOption(silo._DBOPT_TENSOR_RANK, silo._DB_VARTYPE_VECTOR) == 0
    assert optlistVar.addOption(silo._DBOPT_HIDE_FROM_GUI, 1) == 0
    assert optlistVar.addOption(silo._DBOPT_TENSOR_RANK, silo._DB_VARTYPE_SCALAR) == 0

    if dim == 2:
        return ("{%s_x, %s_y}" % (name, name), silo._DB_VARTYPE_VECTOR,
                optlistDef, optlistMV, optlistVar)
    else:
        return ("{%s_x, %s_y, %s_z}" % (name, name, name), silo._DB_VARTYPE_VECTOR,
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
        assert optlist.addOption(silo._DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo._DBOPT_DTIME, time) == 0
    assert optlistMV.addOption(silo._DBOPT_TENSOR_RANK, silo._DB_VARTYPE_TENSOR) == 0
    assert optlistVar.addOption(silo._DBOPT_HIDE_FROM_GUI, 1) == 0
    assert optlistVar.addOption(silo._DBOPT_TENSOR_RANK, silo._DB_VARTYPE_SCALAR) == 0

    if dim == 2:
        return ("{{%s_xx, %s_xy}, {%s_yx, %s_yy}}" % (name, name, name, name), silo._DB_VARTYPE_TENSOR,
                optlistDef, optlistMV, optlistVar)
    else:
        return ("{{%s_xx, %s_xy, %s_xz}, {%s_yx, %s_yy, %s_yz}, {%s_zx, %s_zy, %s_zz}}" % (name, name, name,
                                                                                           name, name, name,
                                                                                           name, name, name),
                silo._DB_VARTYPE_TENSOR, optlistDef, optlistMV, optlistVar)

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
