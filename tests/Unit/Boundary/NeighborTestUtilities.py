from math import *

################################################################################
def distance(point1, point2):
    result = (point1.x - point2.x)**2
    try:
        result = result + (point1.y - point2.y)**2
    except:
        pass
    try:
        result = result + (point1.z - point2.z)**2
    except:
        pass
    result = sqrt(result)
    return result

################################################################################
def findNeighborNodes(position0, radius, nodes):
    result = []
    iNode = 0
    positionField = nodes.positions
    for position in positionField:
        r = distance(position0, position)
        if r <= radius:
            result.append(iNode)
        iNode = iNode + 1
    return result

################################################################################
def checkNeighbors(neighborList, neighborList0):
    if len(neighborList) < len(neighborList0):
        return 0

    result = 1
    i = 0
    while (result == 1 and i < len(neighborList0)):
        if neighborList[:].count(neighborList0[i]) != 1:
            result = 0
        i = i + 1
    return result

################################################################################
def plotNodes(nodes, winNum=0, title=None):
    try:
        from gist import *
        from SpheralGistUtilities import *
        window(winNum)
        x, y = nodePositions2d(nodes)
        plg(y[:nodes.numInternalNodes], x[:nodes.numInternalNodes],
            type=0, marker='\4')
        plg(y[nodes.numInternalNodes:], x[nodes.numInternalNodes:],
            type=0, marker='\4', color='yellow')
        if title:
            pltitle(title)
    except:
        pass
    return

