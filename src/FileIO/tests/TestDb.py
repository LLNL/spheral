from Spheral import *
from SpheralTestUtilities import *

################################################################################
fileName = "test.h5"

file = DbFileIO2d(fileName)
output('file')

################################################################################
title('2-D DB I/O Test: write atomic data types')

groupName = "/Root/AtomicTypes/"

fileInt = 42
fileScalar = 42.5
fileVector = Vector2d(1.3, 4.92)
fileTensor = Tensor2d(0, 1, 2, 3)
fileSymTensor = SymTensor2d(0, 1, 1, 2)
fileString = "Booga booga"

output("fileInt")
file.writeInt(fileInt, groupName + "fileInt")
output("fileScalar")
file.writeScalar(fileScalar, groupName + "fileScalar")
output("fileVector")
file.writeVector(fileVector, groupName + "fileVector")
output("fileTensor")
file.writeTensor(fileTensor, groupName + "fileTensor")
output("fileSymTensor")
file.writeSymTensor(fileSymTensor, groupName + "fileSymTensor")
output("fileString")
file.writeString(fileString, groupName + "fileString")

################################################################################
title('2-D DB I/O Test: write field data types')

groupName = "/Root/FieldTypes/"

n = 50
nodes = SphNodeList2d(n)

scalarField = ScalarField2d(nodes)
vectorField = VectorField2d(nodes)
tensorField = TensorField2d(nodes)
symTensorField = SymTensorField2d(nodes)

for i in xrange(n):
    scalarField[i] = i
    vectorField[i] = (i, 2*i)
    tensorField[i] = Tensor2d(i, 2*i, 3*i, 4*i)
    symTensorField[i] = SymTensor2d(i, 2*i, 2*i, 3*i)

output("scalarField")
file.writeScalarField(scalarField, groupName + "scalarField")
output("vectorField")
file.writeVectorField(vectorField, groupName + "vectorField")
output("tensorField")
file.writeTensorField(tensorField, groupName + "tensorField")
output("symTensorField")
file.writeSymTensorField(symTensorField, groupName + "symTensorField")

################################################################################
title('2-D DB I/O Test: read atomic data types')

# The file is already open for read-write: we don't have to re-open it.

groupName = "/Root/AtomicTypes/"

readInt = file.readInt(groupName + "fileInt")
output("fileInt, readInt")
##assert(fileInt == readInt)

readScalar = file.readScalar(groupName + "fileScalar")
output("fileScalar, readScalar")

readVector = file.readVector(groupName + "fileVector")
output("fileVector, readVector")

readTensor = file.readTensor(groupName + "fileTensor")
output("fileTensor, readTensor")

readSymTensor = file.readSymTensor(groupName + "fileSymTensor")
output("fileSymTensor, readSymTensor")

################################################################################
title('2-D DB I/O Test: read field data types')
################################################################################

groupName = "/Root/FieldTypes/"

readScalarField = ScalarField2d(nodes)
file.readScalarField(readScalarField, groupName + "scalarField")
output("(scalarField[:]), (readScalarField[:])")

readVectorField = VectorField2d(nodes)
file.readVectorField(readVectorField, groupName + "vectorField")
output("(vectorField[:]), (readVectorField[:])")

readTensorField = TensorField2d(nodes)
file.readTensorField(readTensorField, groupName + "tensorField")
output("(tensorField[:]), (readTensorField[:])")

readSymTensorField = SymTensorField2d(nodes)
file.readSymTensorField(readSymTensorField, groupName + "symTensorField")
output("(symTensorField[:]), (readSymTensorField[:])")


