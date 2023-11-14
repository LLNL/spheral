from spheralDimensions import spheralDimensions
dims = spheralDimensions()

def PorousEquationOfState(eosS, *args):
    print("DEPRECATION WARNING: PorousEquationOfState is no longer used, returning unmodified solid equation of state")
    return eosS

for ndim in dims:
    exec(f"PorousEquationOfState{ndim}d = PorousEquationOfState")
