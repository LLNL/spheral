#-------------------------------------------------------------------------------
# Method to find and parse the named enum from a file.
#-------------------------------------------------------------------------------
def parseEnumDef(enumName, fileName):
    f = open(fileName, "r")
    lines = f.readlines()
    result = []
    accumulateLines = False
    for line in lines:
        stuff = line.split()
        if len(stuff) > 0 and accumulateLines and stuff[0] == "};":
            accumulateLines = False
        if len(stuff) > 0 and accumulateLines:
            result.append(stuff[0])
        if len(stuff) > 1 and stuff[0] == "enum" and stuff[1] == enumName:
            accumulateLines = True
    return result
            
#-------------------------------------------------------------------------------
# Add an enum of the given name to the specified scope.
#-------------------------------------------------------------------------------
def addEnumDefinition(scope, enumName, fileName):
    enumValues = parseEnumDef(enumName, fileName)
    return scope.add_enum(enumName, enumValues)
