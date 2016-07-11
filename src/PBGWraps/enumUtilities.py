#-------------------------------------------------------------------------------
# Method to find and parse the named enum from a file.
#-------------------------------------------------------------------------------
def parseEnumDef(enumName, fileName, enclosing_type="enum", val_index=0):
    f = open(fileName, "r")
    lines = f.readlines()
    result = []
    accumulateLines = False
    for line in lines:
        stuff = line.split()
        if len(stuff) > 0 and accumulateLines and stuff[0] == "};":
            accumulateLines = False
        if len(stuff) > 0 and accumulateLines:
            result.append(stuff[val_index])
        if len(stuff) > 1 and stuff[0] == enclosing_type and stuff[1] == enumName:
            accumulateLines = True
    return result
            
#-------------------------------------------------------------------------------
# Add an enum of the given name to the specified scope.
#-------------------------------------------------------------------------------
def addEnumDefinition(scope, enumName, fileName):
    enumValues = parseEnumDef(enumName, fileName)
    return scope.add_enum(enumName, enumValues)

#-------------------------------------------------------------------------------
# Add a struct of the given name with a bunch of integer members.
#-------------------------------------------------------------------------------
def addStructAsEnumDefinition(scope, structName, fileName):
    enumValues = parseEnumDef(structName, fileName, "struct", 3)
    x = scope.add_struct(structName)
    for val in enumValues:
        x.add_static_attribute(val, "long", is_const=True)
    return x
