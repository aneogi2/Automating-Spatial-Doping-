Important Notes:
This function directly modifies the passed DataFrames (DF_atoms, DF_bond, DF_angle). If you need to keep the original DataFrames unchanged, consider making copies before passing them to the function.
The function uses fixed atom types for Na and Cl (4 and 5, respectively). You should adjust these values according to your system's specific atom type assignments.
Hydrogen atoms associated with replaced oxygen atoms are identified by assuming the hydrogen atom IDs are immediately after their corresponding oxygen atom IDs (atom_id + 1). This assumption might not hold in all systems, so adjust the logic if your data structure differs.
The random_state=1 parameter is used for reproducibility in sampling. You can remove or change this value for different random samples.
