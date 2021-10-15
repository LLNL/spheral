from Spheral import *
from SpheralTestUtilities import *

################################################################################
fileName = "test.ascii"
access = 0           # Create
format = 0           # Ascii

file = FlatFileIO2d(fileName, access, format)
output('file')

################################################################################
title('2-D FlatFile I/O Test: write atomic data types')

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
title('2-D FlatFile I/O Test: write Field data types')

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
title('2-D FlatFile I/O Test: read atomic data types')

# Close the current File
file = None

# Open it up again in read only mode.
file = FlatFileIO2d(fileName, 1, format)

groupName = "/Root/AtomicTypes/"

readInt = 0
file.readInt(readInt, groupName + "fileInt")
output("fileInt, readInt")
##assert(fileInt == readInt)

readScalar = 0.0
file.readScalar(readScalar, groupName + "fileScalar")
output("fileScalar, readScalar")

##readVector = Vector2d(0.0, 0.0)
##file.readScalar(readVector, groupName + "fileVector")
##output("fileVector, readVector")

####readTensor = Tensor2d(0, 0, 0, 0)
####file.readTensor(readTensor, groupName + "fileTensor")
####output("fileTensor, readTensor")

####readSymTensor = SymTensor2d(0, 0, 0, 0)
####file.readSymTensor(readSymTensor, groupName + "fileSymTensor")
####output("fileSymTensor, readSymTensor")

readString = ""
file.readString(readString, groupName + "fileString")
output("fileString, readString")

groupName = "/Root/FieldTypes/"

readScalarField = ScalarField2d(nodes)
file.readScalarField(readScalarField, groupName + "scalarField")
output("(scalarField[:]), (readScalarField[:])")

################################################################################
