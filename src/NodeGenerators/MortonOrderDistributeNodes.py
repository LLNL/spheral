from DistributeNodes import makeDistributeNodesMethod
from spheralDimensions import spheralDimensions
dims = spheralDimensions()

for dim in dims:
    exec("import Spheral{dim}d; distributeNodes{dim}d = makeDistributeNodesMethod('{distributor}', Spheral{dim}d)".format(dim=dim, distributor="MortonOrderRedistributeNodes"))
