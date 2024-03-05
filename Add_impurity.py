# Import Modules

import pandas as pd
import numpy as np
import os

'''
  Read the initial structure and disintegrate into Angles Bonds and Positions
'''
def read_angle_data(path_to_file):
    angle_data = []
    filename = f'{path_to_file}/angle_bond'
    with open(filename, "r") as file:
        angle = file.readlines()
    for line in angle:
        ang_col = line.split()
        if len(ang_col) == 5:
            ang_col = [int(val) for val in ang_col]
            angle_data.append(ang_col)
    DF_angle = pd.DataFrame(np.array(angle_data), columns=["angle_id", "angle_type", "H", "O", "H"])
    DF_angle.to_csv("angles.csv", sep=' ', index=False)
    return DF_angle

def write_angles(Natoms):
    Natoms = 2880  # Assuming this is defined or passed as an argument
    Nbonds = int((2*Natoms)/3)
    Nangles = int((Nbonds/2))
    initial_angle_data = {'angle_id': [1], 'angle_type': [1], '1st_H': [2], 'O': [1], '2nd_H': [3]}
    df_initial = pd.DataFrame(initial_angle_data)
    df_initial.to_csv('init_angle.csv', sep=' ', index=False)

#Write the Bonds and atoms connecting to the bonds based on the Angle information

def write_bonds(Natoms, filename_output):
    Nbonds = int((2*Natoms)/3)

    Bdf = pd.DataFrame(index=range(Nbonds+1), columns=["bond_id", "bond_type", "O", "H"])
    Bdf["bond_id"] = Bdf.index
    Bdf["bond_type"] = 1
    Bdf = Bdf.drop([0])

    Bdf.iloc[0, 2] = 1  # 1st Oxygen atom ID
    Bdf.iloc[0, 3] = 2  # 1st Hydrogen atom ID
    Bdf.iloc[1, 2] = 1
    Bdf.iloc[1, 3] = 3

    for i in range(3, len(Bdf)):
        if Bdf.loc[i-1, 'O'] == Bdf.loc[i-2, 'O']:
            Bdf.loc[i, 'O'] = Bdf.loc[i-1, 'O'] + 3
        else:
            Bdf.loc[i, 'O'] = Bdf.loc[i-1, 'O']
    Bdf['O'].fillna(method='ffill', inplace=True)

    for i in range(3, len(Bdf)):
        if Bdf.loc[i-1, 'H'] == Bdf.loc[i-2, 'H'] + 1:
            Bdf.loc[i, 'H'] = Bdf.loc[i-1, 'H'] + 2
        else:
            Bdf.loc[i, 'H'] = Bdf.loc[i-1, 'H'] + 1
    Bdf['H'].iloc[-1] = int(Natoms)

    Bdf.to_csv(filename_output, sep=' ', index=False, header=None)

# Reset atom ids based combining the bonds and  angles

def reset_atom_ids(filename_input, filename_output):
    data = []
    with open(filename_input, "r") as file:
        atoms = file.readlines()
    for line in atoms:
        atom_col = line.split()
        if len(atom_col) > 6:
            atom_col = [int(val) if i != 3 else float(val) for i, val in enumerate(atom_col)]
            data.append(atom_col)

    df_data = pd.DataFrame(data, columns=["atom_id", "mol_id", "atom_type", "charge", "x", "y", "z"])
    df_data = df_data.drop([0])
    df_data["atom_id"] = df_data.index

    df_data.to_csv(filename_output, sep=' ', index=False)

# Combine the Bonds Angles and positional information for the pure Ice structure and the read the pure ice structure separately

def read_ice_data(base_path):
    # Initialize empty lists for data
    ice_data = []
    bond_data = []
    angle_data = []

    # Path to the ice data file
    filename = f'{base_path}/with_tip4p_cutoff/with_cu/Ih_cu.geo'
    
    # Read atom data
    with open(filename, "r") as file:
        for line in file:
            atom_col = line.split()
            if len(atom_col) == 7:
                atom_col = [float(val) for val in atom_col]
                ice_data.append(atom_col)
    
    # Convert atom data to DataFrame
    DF_atoms = pd.DataFrame(ice_data, columns=["atom_id", "mol_id", "atom_type", "charge", "x", "y", "z"])
    
    # Read bond data
    with open(filename, "r") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if "Bonds" in line:
                for bline in lines[i+2:]:
                    bcol = bline.split()
                    if len(bcol) == 4:
                        bcol = [int(val) for val in bcol]
                        bond_data.append(bcol)
                break
    
    # Convert bond data to DataFrame
    DF_bond = pd.DataFrame(bond_data, columns=["Bond_id", "Bond_type", "O_id", "H_id"])
    
    # Read angle data
    with open(filename, "r") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if "Angles" in line:
                for aline in lines[i+2:]:
                    acol = aline.split()
                    if len(acol) == 5:
                        acol = [int(val) for val in acol]
                        angle_data.append(acol)
                break
    
    # Convert angle data to DataFrame
    DF_angle = pd.DataFrame(angle_data, columns=["Angle_id", "Angle_type", "H_id", "O_id", "H2_id"])
    
    return DF_atoms, DF_bond, DF_angle
'''
Add any impururites/dopants to the Ice structure in any specific spatial region
'''
# Adding "SALT TO INJURY"

def add_salt_to_ice(DF_atoms, DF_bond, DF_angle, num_of_atoms, Wt_pct_salt=1):
    # Define atomic masses
    mass_H = 1.008
    mass_O = 15.9994
    mass_Na = 22.989
    mass_Cl = 35.453

    # Calculate total initial mass of the system
    init_mass = num_of_atoms * (mass_H + mass_O)

    # Calculate mass and number of NaCl units required
    NaCl_mass = (Wt_pct_salt / 100) * init_mass
    num_NaCl_atoms = round(NaCl_mass / (mass_Na + mass_Cl))
    num_Na_atoms = num_NaCl_atoms // 2
    num_Cl_atoms = num_NaCl_atoms - num_Na_atoms

    # Replace oxygen atoms with Na and Cl
    DF_Na = DF_atoms.query('z < 18 & z > 10 & atom_type == 2.0').sample(num_Na_atoms, random_state=1)
    DF_Na['atom_type'] = 4  # Assuming atom type 4 is Na
    DF_Cl = DF_atoms.query('z < 18 & z > 10 & atom_type == 2.0').sample(num_Cl_atoms, random_state=1)
    DF_Cl['atom_type'] = 5  # Assuming atom type 5 is Cl

    # Update DF_atoms with Na and Cl
    DF_atoms.update(DF_Na)
    DF_atoms.update(DF_Cl)

    # Identify hydrogen atoms to delete (assuming H atoms have atom_type == 1)
    H_to_delete = DF_atoms[DF_atoms['atom_type'].isin([4, 5])]['atom_id'].values + 1
    DF_atoms = DF_atoms[~DF_atoms['atom_id'].isin(H_to_delete)]

    # Update bond data
    DF_bond = DF_bond[~DF_bond['O_id'].isin(DF_Na['atom_id']) & ~DF_bond['O_id'].isin(DF_Cl['atom_id'])]

    # Update angle data
    DF_angle = DF_angle[~DF_angle['O_id'].isin(DF_Na['atom_id']) & ~DF_angle['O_id'].isin(DF_Cl['atom_id'])]

    return DF_atoms, DF_bond, DF_angle


