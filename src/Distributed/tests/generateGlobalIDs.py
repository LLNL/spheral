#===============================================================================
# Generate unique global IDs for a set of NodeLists.
#===============================================================================
def generateGlobalIDs(nodeLists,
                      globalNodeIDs,
                      numGlobalNodes):
    result = []
    offset = 0
    for nodes in nodeLists:
        ids = globalNodeIDs(nodes)
        for i in xrange(nodes.numInternalNodes):
            ids[i] += offset
        offset += numGlobalNodes(nodes)
        result.append(ids)
    return tuple(result)

