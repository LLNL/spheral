from spheralDimensions import spheralDimensions
dims = spheralDimensions()

def PorousStrengthModel(strengthModelS, *args):
    print("DEPRECATION WARNING: PorousStrengthModel is no longer used")
    return strengthModelS

for ndim in dims:
    exec(f"PorousStrengthModel{ndim}d = PorousStrengthModel")
