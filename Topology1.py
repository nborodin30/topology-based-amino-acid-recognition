from itertools import islice
import sys

# Get file path from command-line argument
if len(sys.argv) < 2:
    print("Usage: python Topology1.py <sdf_file_path>")
    sys.exit(1)
file_path = sys.argv[1]

# Read all lines of the raw data
with open(file_path, "r") as sdf_file:
    raw_data = sdf_file.readlines()

'''
Extracts atoms and bonds from the SDF file.
Returns a dictionary with atoms as keys and their bonds as values.
'''
def get_atoms_and_bond(df_sdf):
    atoms = []
    atom_section = False
    bonds = []
    atoms_with_bonds = {}
    bond_section = False
    
    # Parse the SDF file to extract atoms and bonds
    for line in raw_data:
        if line.startswith("M  V30 BEGIN ATOM"):
            atom_section = True
            continue
        if line.startswith("M  V30 END ATOM"):
            atom_section = False
            continue
        if line.startswith("M  V30 BEGIN BOND"):
            bond_section = True
            continue
        if line.startswith("M  V30 END BOND"):
            bond_section = False
            continue
        
        # Extract atom information
        if atom_section and line.startswith("M  V30 "):
            atom = line.split()[3]
            atoms.append(atom)
            
        # Extract bond information
        if bond_section and line.startswith("M  V30 "):
            bond_data = line.split()
            id = int(bond_data[2])
            type_of_bond = int(bond_data[3])
            atom1 = int(bond_data[4])  
            atom2 = int(bond_data[5])
            bonds.append((id,type_of_bond,atom1, atom2))
            
    # Create a dictionary linking atoms to their bonds
    for i in range(len(atoms)):
        atoms_with_bonds[atoms[i] + " "+ str(i+1)] = []
        for bond in bonds:
            if i+1 == bond[2]:
                atoms_with_bonds[atoms[i] + " "+ str(i+1)].append(atoms[bond[3]-1] + " id: " + str(bond[3]) + " Type of bond: " + str(bond[1]))
            if i+1 == bond[3]:
                atoms_with_bonds[atoms[i] + " "+ str(i+1)].append(atoms[bond[2]-1] + " id: "+ str(bond[2])+ " Type of bond: " + str(bond[1]))
                
    # Remove atoms with no bonds
    for key in list(atoms_with_bonds.keys()):
        if atoms_with_bonds[key] == []:
            del atoms_with_bonds[key]
        
    return atoms_with_bonds

'''
Identifies NCCO patterns (peptide backbone) for a given atom.
Returns a list of atoms forming the NCCO pattern if found, or None.
'''
def check_NCCO(atom, atoms_with_bonds): 
    atom_and_id = atom.split()
    result = []
    NC= CC = CO = False
    
    # Check if the atom is Nitrogen and has a bond to a Carbon
    if atom_and_id[0] == "N":
        for bond in atoms_with_bonds[atom]:
            bond_and_id = bond.split() 
            if bond_and_id[0] == "C" and bond_and_id[6] == "1" and int(bond_and_id[2]) == int(atom_and_id[1])+1:
                NC = True
                result.append(atom)
                break
                
    # If N-C bond found, check for C-C bond
    if NC:
        for bond in atoms_with_bonds[f'C {int(atom_and_id[1]) +1 }']:
            bond_and_id = bond.split()               
            if bond_and_id[0] == "C" and bond_and_id[6] == "1" and int(bond_and_id[2]) == int(atom_and_id[1])+2:
                CC = True
                result.append(f'C {int(atom_and_id[1])+1}')
                break
                
    # If C-C bond found, check for C=O bond
    if CC:
        for bond in atoms_with_bonds[f'C {int(atom_and_id[1])+2}']:
            bond_and_id = bond.split()               
            if bond_and_id[0] == "O" and bond_and_id[6] == "2" and int(bond_and_id[2]) == int(atom_and_id[1])+3:
                CO = True
                result.append(f'C {int(atom_and_id[1])+2}')
                result.append(f'O {int(atom_and_id[1])+3}')
                break
                
    # Return the NCCO pattern if found
    if NC and CC and CO:
        return result    

'''
Finds all NCCO patterns (peptide backbones) in the molecule.
Returns a list of all atoms that are part of a peptide backbone.
'''
def all_NCCO(atoms_with_bonds):
    nccos = []
    for atom in atoms_with_bonds.keys():
        ncco = check_NCCO(atom, atoms_with_bonds)
        if ncco:
            nccos.extend(ncco)
    return nccos

'''
Finds the last oxygen atom in the backbone with a single bond.
Used to determine where the peptide chain ends.
'''
def last_atom(atom_bonds, backbone):
    for atom in atom_bonds.keys():
        if atom.startswith('O'):
            for bond in atom_bonds[atom]:
                if bond.split()[0] + " " + bond.split()[2] in backbone and bond.split()[6] == '1':
                    return atom

'''
Separates atoms into groups corresponding to amino acids.
Returns a dictionary where each key is a group ID and each value is a dictionary
of atoms and their bonds in that amino acid.
'''
def side_chains(atoms_with_bonds, backbone):
    groups = {}        # Dictionary for groups
    rest_for_ac = {}   # Dictionary for the current group
    current_group = [] # Tracks atoms in the current NCCO group
    counter = 0        # Group counter
    
    # Get all backbone atoms and the last atom of the chain
    ncco_atoms = set(backbone)
    last_backbone_atom = last_atom(atoms_with_bonds, backbone)
    
    for atom in atoms_with_bonds.keys():
        # If atom is part of the backbone
        if atom in ncco_atoms:
            # If atom is an oxygen, save current side chain group and reset
            if atom.startswith('O'):
                groups[f'group_{counter}'] = rest_for_ac
                rest_for_ac = {}  # Reset for the next group
                
            # Add atom to current NCCO group
            current_group.append(atom)
            
            # If we've completed a full NCCO pattern, increment counter
            if len(current_group) == 4:
                counter += 1
                current_group = []  # Reset the current NCCO group
        else:
            # If atom is not in backbone, add to current side chain
            rest_for_ac[atom] = atoms_with_bonds[atom]
            
        # If we've reached the end of the chain, save the last group and return
        if atom == last_backbone_atom:
            rest_for_ac.pop(last_backbone_atom)
            groups[f'group_{counter}'] = rest_for_ac
            return groups
            
    # Save the last group and return
    groups[f'group_{counter}'] = rest_for_ac
    return groups

'''
Translates long amino acid names to their three-letter codes.
'''
def translate_AC_to_name(AC):
    bib = {
        'Alanine': 'ALA',
        'Arginine': 'ARG',
        'Asparagine': 'ASN',
        'Aspartate': 'ASP',
        'Cysteine': 'CYS',
        'Glutamate': 'GLU',
        'Glutamine': 'GLN',
        'Glycine': 'GLY',
        'Histidine': 'HIS',
        'Isoleucine': 'ILE',
        'Leucine': 'LEU',
        'Lysine': 'LYS',
        'Methionine': 'MET',
        'Phenylalanine': 'PHE',
        'Proline': 'PRO',
        'Serine': 'SER',
        'Threonine': 'THR',
        'Tryptophan': 'TRP',
        'Tyrosine': 'TYR',
        'Valine': 'VAL',
        'Unknown' : 'UNK',
        }
    return bib[AC]

'''
Identifies the amino acid type based on the side chain structure.
Returns the name of the identified amino acid.
'''
def side_chain_To_AC(group, backbone):
    AC= "Unknown"
    first_atom_in_chain = find_first_after_C_alpha(group, backbone)
    
    # If no side chain atoms found, it's Glycine
    if not first_atom_in_chain:
        AC = "Glycine"
        return AC
        
    # Analyze bonds of the first atom in the side chain
    first_score = score_for_bonds(first_atom_in_chain, group)
    
    # Identify based on first atom's bonds
    if first_score['C'][0] == 1 and first_score['C'][1] == '1' and len(group) == 1:
        AC = "Alanine"
        return AC
    if first_score['C'][0] == 1 and first_score['C'][1] == '1' and first_score['O'][0] == 1 and first_score['O'][1] == '1':
        AC = "Serine"
        return AC
    if first_score['C'][0] == 2 and first_score['C'][1] == '11' and first_score['O'][0] == 1 and first_score['O'][1] == '1':
        AC = "Threonine"
        return AC
    if first_score['C'][0] == 1 and first_score['C'][1] == '1' and first_score['S'][0] == 1 and first_score['S'][1] == '1':
        AC = "Cysteine"
        return AC
        
    # Move to second atom if present
    second_atom_in_chain = one_step(first_atom_in_chain, first_atom_in_chain, group, backbone)
    if second_atom_in_chain != None:
        second_score = score_for_bonds(second_atom_in_chain, group)
        
        # Identify based on second atom's bonds
        if second_score['C'][0] == 3 and second_score['C'][1] == '111':
            AC = "Leucine"
            return AC
        if second_score['C'][0] == 2 and second_score['C'][1] == '11' and len(group) == 4:
            AC = "Isoleucine"
            return AC
        if second_score['S'][0] == 1 and second_score['S'][1] == '1':
            AC = "Methionine"
            return AC
        if second_score['C'][0] == 1 and second_score['C'][1] == '1' and second_score['O'][0] == 1 and second_score['O'][1] == '2' and second_score['N'][0] == 1 and second_score['N'][1] == '1':   
            AC = "Asparagine"
            return AC
        if second_score['C'][0] == 1 and second_score['C'][1] == '1' and second_score['O'][0] == 2 and (second_score['O'][1] == '12' or second_score['O'][1] == '21'):   
            AC = "Aspartate"
            return AC
        if second_score['C'][0] == 2 and (second_score['C'][1] == '12' or second_score['C'][1] =='21' or second_score['C'][1]=='11') and second_score['N'][0] == 1 and (second_score['N'][1] == '1' or second_score['N'][1]=='2'):
            AC = "Histidine"
            return AC
        if second_score['C'][0] == 3 and (second_score['C'][1] == '112' or second_score['C'][1] == '121' or second_score['C'][1] =='211'):
            if len(group) == 7:
                AC = "Phenylalanine"
                return AC
            if len(group) == 8:
                AC = "Tyrosine"
                return AC
            if len(group) == 10:
                AC = "Tryptophan"
                return AC
    else: 
        AC = "Valine"
        return AC
        
    # Move to third atom if present
    third_atom_in_chain = one_step(first_atom_in_chain, second_atom_in_chain, group, backbone)
    if third_atom_in_chain != None:
        third_score = score_for_bonds(third_atom_in_chain, group)
        
        # Identify based on third atom's bonds
        if third_score['C'][0] == 1 and third_score['C'][1] == '1' and third_score['O'][0] == 2 and (third_score['O'][1] == '21' or third_score['O'][1] =='12'):
            AC = "Glutamate"
            return AC
        if third_score['C'][0] == 1 and (third_score['C'][1] == '1' or third_score['C'][1] == '2') and third_score['N'][0] == 1 and third_score['N'][1] == '1':
            for bond in group[third_atom_in_chain]:
                if bond.split()[0] == 'N' and bond.split()[6] == '1' and ((bond.split()[0] + " " + bond.split()[2]) in backbone):
                    AC = "Proline"
                    return AC
        if third_score['C'][0] == 1 and third_score['C'][1] == '1' and third_score['O'][0] == 1 and third_score['O'][1] == '2' and third_score['N'][0] == 1 and third_score['N'][1] == '1':
            AC = "Glutamine"
            return AC
            
    # Move to fourth atom if present
    fourth_atom_in_chain = one_step(second_atom_in_chain, third_atom_in_chain, group, backbone)
    if fourth_atom_in_chain != None:
        fourth_score = score_for_bonds(fourth_atom_in_chain, group)
        
        # Identify based on fourth atom's bonds
        if fourth_score['C'][0] == 1 and fourth_score['C'][1] == '1' and fourth_score['N'][0] == 1 and fourth_score['N'][1] == '1':
            AC = "Lysine"
            return AC
        if fourth_score['C'][0] == 2 and (fourth_score['C'][1] == '11' or fourth_score['C'][1] == '12' or fourth_score['C'][1] == '21'):
            AC = "Arginine"
            return AC
            
    return AC

'''
Finds the first atom after C-alpha in the side chain.
Returns the key for this atom.
'''
def find_first_after_C_alpha(bonds, backbone):
    C_alpha = None
    for atom in bonds.keys():
        if atom.startswith('C'):
            for bond in bonds[atom]:
                C_alpha_name = bond.split()[0] + " " + bond.split()[2]
                if C_alpha_name in backbone and C_alpha_name.startswith('C'):
                    C_alpha = atom
                    break
    return C_alpha

'''
Calculates the number and types of bonds for a given atom.
Returns a dictionary with counts and bond types.
'''
def score_for_bonds(atom, group):
    score = {'C': [0,''], 'O': [0,''], 'N': [0,''] ,'S': [0,'']}
    for bond in group[atom]:
        atom_type = bond.split()[0]
        atom_bond_type = bond.split()[6]
        if atom_type == 'C':
            score['C'][0] += 1
            score['C'][1] += atom_bond_type
        if atom_type == 'O':
            score['O'][0] += 1
            score['O'][1] += atom_bond_type
        if atom_type == 'N':
            score['N'][0] += 1
            score['N'][1] += atom_bond_type
        if atom_type == 'S':
            score['S'][0] += 1
            score['S'][1] += atom_bond_type
    return score

'''
Finds the next atom in a chain from a previous atom.
Used for traversing side chains to identify amino acids.
'''
def one_step(first_atom, prev_atom, group, backbone):
    for atom in group.keys():
        if atom not in backbone:
            for bond in group[atom]:
                if bond.split()[0] + " " + bond.split()[2] == prev_atom and len(group[atom]) > 1 and atom != first_atom:
                    return atom
    return None

'''
Maps each atom to its corresponding amino acid group.
Returns a dictionary mapping atom IDs to group IDs.
'''
def assign_atoms_to_groups(sides, backbone, atoms_with_bonds):
    atom_to_group = {}
    
    # First assign side chain atoms to their groups
    for group_id, group in sides.items():
        for atom_key in group.keys():
            atom_id = int(atom_key.split()[1])
            atom_to_group[atom_id] = group_id
    
    # Then assign backbone atoms to groups
    ncco_groups = {}
    current_group = None
    for atom in backbone:
        atom_id = int(atom.split()[1])
        if atom.startswith('N'):  # Start of new NCCO group
            current_group = (atom_id // 4) + 1  # Approximate grouping
        
        if current_group is not None:
            atom_to_group[atom_id] = f'group_{current_group}'
    
    return atom_to_group

'''
Creates an annotated SDF file with amino acid labels for each atom.
Writes the output to the specified file.
'''
def create_annotated_sdf(input_data, output_file, sides, backbone, atoms_with_bonds):
    # Map amino acids to groups directly
    group_to_aa = {}
    for group_id, group in sides.items():
        group_to_aa[group_id] = translate_AC_to_name(side_chain_To_AC(group, backbone))
    
    # Map atoms to groups
    atom_to_group = assign_atoms_to_groups(sides, backbone, atoms_with_bonds)
    
    # Create the annotated SDF file
    with open(output_file, "w") as outfile:
        in_atom_section = False
        
        for line in input_data:
            # Check if we're in the atom section
            if line.startswith("M  V30 BEGIN ATOM"):
                in_atom_section = True
                outfile.write(line)
                continue
            
            if line.startswith("M  V30 END ATOM"):
                in_atom_section = False
                outfile.write(line)
                continue
            
            # If in atom section, append amino acid
            if in_atom_section and line.startswith("M  V30") and not line.startswith("M  V30 BEGIN") and not line.startswith("M  V30 END"):
                parts = line.split()
                if len(parts) >= 3:
                    atom_id = int(parts[2])
                    
                    # Find amino acid for this atom
                    aa_code = "UNK"  # Default
                    if atom_id in atom_to_group:
                        group_id = atom_to_group[atom_id]
                        if group_id in group_to_aa:
                            aa_code = group_to_aa[group_id]
                    
                    # Append amino acid to line
                    outfile.write(f"{line.rstrip()} {aa_code}\n")
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

# Main execution
atoms_with_bonds = get_atoms_and_bond(raw_data)
backbone = all_NCCO(atoms_with_bonds)
sides = side_chains(atoms_with_bonds, backbone)
if 'group_0' in sides:
    sides.pop('group_0')

# Create annotated SDF file
import os
base_name = os.path.splitext(os.path.basename(file_path))[0]
output_sdf_path = f"{base_name}-annotated.sdf"
print(f"\nCreating annotated SDF file at {output_sdf_path}")
create_annotated_sdf(raw_data, output_sdf_path, sides, backbone, atoms_with_bonds)
