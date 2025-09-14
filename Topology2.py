#import all modules that would be used
from itertools import islice
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from grakel import Graph, GraphKernel

import sys
# Get file path from command-line argument
if len(sys.argv) < 2:
    print("Usage: python Topology2.py <sdf_file_path>")
    sys.exit(1)
file_path = sys.argv[1]

# Read all lines of the raw data
with open(file_path, "r") as sdf_file:
    raw_data = sdf_file.readlines()

'''
Parses the SDF data to create a NetworkX graph with atoms as nodes and bonds as edges.
Bond weights represent bond types: single_bond -> '1', double_bond -> '2', triple_bond -> '3'.
'''
def get_nodes_and_edges(data):
    G = nx.Graph()
    atom_section = False
    bond_section = False
    for line in data:
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
        if atom_section:
            atom = line.split()
            G.add_node(atom[2], atom=atom[3])
        if bond_section:
            bond = line.split()
            node1 = bond[4]
            node2 = bond[5]
            G.add_edge(node1, node2, bond_type=bond[3])
    return G

'''
Identifies peptide backbone atoms (NCCO pattern) and C-alpha atoms.
Returns two lists:
1. All atoms that are part of the backbone
2. Only the C-alpha atoms
'''
def find_all_backbones_and_C_alpha(G):  
    backbones = []
    c_alphas = []
    for n_node in G.nodes():
        if G.nodes[n_node].get('atom') != 'N':
            continue

        for ca_node in G.neighbors(n_node):
            if G.nodes[ca_node].get('atom') != 'C':
                continue
            edge_data = G.get_edge_data(n_node, ca_node)
            if not (edge_data and edge_data.get('bond_type') == '1'):
                continue
            
            for c_node in G.neighbors(ca_node):
                if c_node == n_node or G.nodes[c_node].get('atom') != 'C':
                    continue
                edge_data_ca_c = G.get_edge_data(ca_node, c_node)
                if not (edge_data_ca_c and edge_data_ca_c.get('bond_type') == '1'):
                    continue

                for o_node in G.neighbors(c_node):
                    if o_node in (ca_node, n_node) or G.nodes[o_node].get('atom') != 'O':
                        continue
                    edge_data_c_o = G.get_edge_data(c_node, o_node)
                    if not (edge_data_c_o and edge_data_c_o.get('bond_type') == '2'):
                        continue

                    backbones.append(n_node)
                    backbones.append(ca_node)
                    c_alphas.append(ca_node)
                    backbones.append(c_node)
                    backbones.append(o_node)   
    return backbones, c_alphas

'''
Extracts the side chain graph for a given C-alpha atom.
Returns a NetworkX graph representing the side chain.
'''
def extract_side_chain(c_alpha, backbone, C_alphas, Graphh):
    Gr = nx.Graph()
    visited = set()

    def recurse(atom):
        if atom in visited:
            return
        visited.add(atom)

        if not Gr.has_node(atom):
            Gr.add_node(atom, atom=Graphh.nodes[atom]['atom'])

        for neighbor in Graphh.neighbors(atom):
            try:
                neighbor_atom_type = Graphh.nodes[neighbor]['atom']
                # Only go deeper if neighbor is NOT in the backbone
                if neighbor not in backbone or (neighbor_atom_type == 'N' and atom not in C_alphas):
                    if not Gr.has_node(neighbor):
                        Gr.add_node(neighbor, atom=neighbor_atom_type)

                    bond_data = Graphh.get_edge_data(atom, neighbor)
                    if not Gr.has_edge(atom, neighbor):
                        Gr.add_edge(atom, neighbor, bond_type=bond_data['bond_type'])

                    recurse(neighbor)
            except Exception as e:
                continue

    recurse(c_alpha)
    return Gr

'''
Creates a graph pattern for alanine side chain.
'''
def alanine():
    side_chain_alanine_graph = nx.Graph()
    side_chain_alanine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_alanine_graph.add_node('C2', atom='C')  # Carbon in side chain (CH3)
    side_chain_alanine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    return side_chain_alanine_graph

'''
Creates a graph pattern for glycine side chain.
'''
def glycine():
    side_chain_glycine_graph = nx.Graph()
    side_chain_glycine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    return side_chain_glycine_graph

'''
Creates a graph pattern for valine side chain.
'''
def valine():
    side_chain_valine_graph = nx.Graph()
    side_chain_valine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_valine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_valine_graph.add_node('C3', atom='C')  # Second carbon in side chain
    side_chain_valine_graph.add_node('C4', atom='C')  # Third carbon in side chain (CH3)
    
    side_chain_valine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond from Cα to C2
    side_chain_valine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond from C2 to C3
    side_chain_valine_graph.add_edge('C2', 'C4', bond_type='1')  # Single bond from C2 to C4 (CH3 group)
    
    return side_chain_valine_graph

'''
Creates a graph pattern for leucine side chain.
'''
def leucine():
    side_chain_leucine_graph = nx.Graph()
    
    # Define nodes (atoms)
    side_chain_leucine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_leucine_graph.add_node('C2', atom='C')  # Carbon in side chain (CH2)
    side_chain_leucine_graph.add_node('C3', atom='C')  # Second carbon in side chain (CH)
    side_chain_leucine_graph.add_node('C4', atom='C')  # Third carbon in side chain (CH2)
    side_chain_leucine_graph.add_node('C5', atom='C')  # Fourth carbon in side chain (CH2)
    
    # Add edges with bond types
    side_chain_leucine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_leucine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_leucine_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond
    side_chain_leucine_graph.add_edge('C3', 'C5', bond_type='1')  # Single bond

    return side_chain_leucine_graph

'''
Creates a graph pattern for isoleucine side chain.
'''
def isoleucine():
    side_chain_isoleucine_graph = nx.Graph()
    side_chain_isoleucine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_isoleucine_graph.add_node('C1', atom='C')  # Carbon in side chain
    side_chain_isoleucine_graph.add_node('C2', atom='C')  # Additional carbon in side chain
    side_chain_isoleucine_graph.add_node('C3', atom='C')  # Additional carbon in side chain
    side_chain_isoleucine_graph.add_node('C4', atom='C')  # Additional carbon in side chain
    
    side_chain_isoleucine_graph.add_edge('Cα', 'C1', bond_type='1')  # Single bond
    side_chain_isoleucine_graph.add_edge('C1', 'C3', bond_type='1')  # Single bond
    side_chain_isoleucine_graph.add_edge('C1', 'C2', bond_type='1')  # Single bond
    side_chain_isoleucine_graph.add_edge('C2', 'C4', bond_type='1')  # Single bond
    return side_chain_isoleucine_graph

'''
Creates a graph pattern for serine side chain.
'''
def serine():
    side_chain_serine_graph = nx.Graph()
    side_chain_serine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_serine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_serine_graph.add_node('O', atom='O')   # Oxygen in side chain
    
    side_chain_serine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_serine_graph.add_edge('C2', 'O', bond_type='1')   # Single bond
    return side_chain_serine_graph

'''
Creates a graph pattern for cysteine side chain.
'''
def cysteine():
    side_chain_cysteine_graph = nx.Graph()
    side_chain_cysteine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_cysteine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_cysteine_graph.add_node('S', atom='S')   # Sulfur in side chain
    side_chain_cysteine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_cysteine_graph.add_edge('C2', 'S', bond_type='1')   # Single bond
    return side_chain_cysteine_graph

'''
Creates a graph pattern for methionine side chain.
'''
def methionine():
    side_chain_methionine_graph = nx.Graph()
    side_chain_methionine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_methionine_graph.add_node('C2', atom='C')  # Carbon in side chain  
    side_chain_methionine_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_methionine_graph.add_node('S', atom='S')  # Sulfur in side chain
    side_chain_methionine_graph.add_node('C4', atom='C')  # Additional carbon in side chain

    side_chain_methionine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_methionine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_methionine_graph.add_edge('C3', 'S', bond_type='1')   # Single bond
    side_chain_methionine_graph.add_edge('S', 'C4', bond_type='1')    # Single bond
    return side_chain_methionine_graph

'''
Creates a graph pattern for threonine side chain.
'''
def threonine():
    side_chain_threonine_graph = nx.Graph()
    side_chain_threonine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_threonine_graph.add_node('C2', atom='C')
    side_chain_threonine_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_threonine_graph.add_node('O', atom='O')   # Hydroxyl group in side chain
    
    side_chain_threonine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_threonine_graph.add_edge('C2', 'O', bond_type='1')   # Single bond
    side_chain_threonine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    return side_chain_threonine_graph

'''
Creates a graph pattern for aspartate side chain.
'''
def aspartate():
    side_chain_aspartate_graph = nx.Graph()
    side_chain_aspartate_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_aspartate_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_aspartate_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_aspartate_graph.add_node('O', atom='O')   # Hydroxyl group
    side_chain_aspartate_graph.add_node('O2', atom='O')  # Second oxygen in side chain
    
    side_chain_aspartate_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_aspartate_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_aspartate_graph.add_edge('C3', 'O', bond_type='2')   # Double bond (carboxyl group)
    side_chain_aspartate_graph.add_edge('C3', 'O2', bond_type='1')  # Single bond (carboxyl group)
    return side_chain_aspartate_graph

'''
Creates a graph pattern for asparagine side chain.
'''
def asparagine():
    side_chain_asparagine_graph = nx.Graph()
    side_chain_asparagine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_asparagine_graph.add_node('C2', atom='C') 
    side_chain_asparagine_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_asparagine_graph.add_node('O', atom='O')   # Oxygen in side chain (double bond in amide group)
    side_chain_asparagine_graph.add_node('N', atom='N')   
    
    side_chain_asparagine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_asparagine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_asparagine_graph.add_edge('C3', 'O', bond_type='2')  # Double bond (amide group)
    side_chain_asparagine_graph.add_edge('C3', 'N', bond_type='1')  # Single bond (amide group)
    return side_chain_asparagine_graph

'''
Creates a graph pattern for glutamate side chain.
'''
def glutamate():
    side_chain_glutamate_graph = nx.Graph()
    side_chain_glutamate_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_glutamate_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_glutamate_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_glutamate_graph.add_node('C4', atom='C')  # Carbon in side chain
    side_chain_glutamate_graph.add_node('O', atom='O')   # Hydroxyl group
    side_chain_glutamate_graph.add_node('O2', atom='O')  # Second oxygen in side chain
    
    side_chain_glutamate_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_glutamate_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_glutamate_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond
    side_chain_glutamate_graph.add_edge('C4', 'O', bond_type='2')   # Double bond (carboxyl group)
    side_chain_glutamate_graph.add_edge('C4', 'O2', bond_type='1')  # Single bond (carboxyl group)
    return side_chain_glutamate_graph

'''
Creates a graph pattern for glutamine side chain.
'''
def glutamine():
    side_chain_glutamine_graph = nx.Graph()
    side_chain_glutamine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_glutamine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_glutamine_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_glutamine_graph.add_node('C4', atom='C')  # Carbon in side chain
    side_chain_glutamine_graph.add_node('O', atom='O')   # Oxygen in side chain (double bond in amide group)
    side_chain_glutamine_graph.add_node('N', atom='N')   # Nitrogen in side chain (amide group)

    side_chain_glutamine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_glutamine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_glutamine_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond
    side_chain_glutamine_graph.add_edge('C4', 'O', bond_type='2')   # Double bond (amide group)
    side_chain_glutamine_graph.add_edge('C4', 'N', bond_type='1')   # Single bond (amide group)
    return side_chain_glutamine_graph

'''
Creates a graph pattern for lysine side chain.
'''
def lysine():
    side_chain_lysine_graph = nx.Graph()
    side_chain_lysine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_lysine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_lysine_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_lysine_graph.add_node('C4', atom='C')  # Carbon in side chain
    side_chain_lysine_graph.add_node('C5', atom='C')  # Carbon in side chain
    side_chain_lysine_graph.add_node('N', atom='N')   # Nitrogen in side chain
    
    side_chain_lysine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_lysine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_lysine_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond
    side_chain_lysine_graph.add_edge('C4', 'C5', bond_type='1')  # Single bond
    side_chain_lysine_graph.add_edge('C5', 'N', bond_type='1')   # Single bond (amine group)
    return side_chain_lysine_graph

'''
Creates a graph pattern for arginine side chain.
'''
def arginine():
    side_chain_arginine_graph = nx.Graph()
    side_chain_arginine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_arginine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_arginine_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_arginine_graph.add_node('C4', atom='C')  # Carbon in side chain
    side_chain_arginine_graph.add_node('N', atom='N')   # Nitrogen in side chain
    side_chain_arginine_graph.add_node('C5', atom='C')  
    side_chain_arginine_graph.add_node('N2', atom='N')  
    side_chain_arginine_graph.add_node('N3', atom='N')  
    
    side_chain_arginine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_arginine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_arginine_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond 
    side_chain_arginine_graph.add_edge('C4', 'N', bond_type='1')   # Single bond 
    side_chain_arginine_graph.add_edge('N', 'C5', bond_type='1')   # Single bond 
    side_chain_arginine_graph.add_edge('C5', 'N2', bond_type='2')  # Double bond 
    side_chain_arginine_graph.add_edge('C5', 'N3', bond_type='1')  # Single bond 
    return side_chain_arginine_graph

'''
Creates another graph pattern for arginine side chain (isomer).
'''
def arginine1():
    side_chain_arginine_graph = nx.Graph()
    side_chain_arginine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_arginine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_arginine_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_arginine_graph.add_node('C4', atom='C')  # Carbon in side chain
    side_chain_arginine_graph.add_node('N', atom='N')   # Nitrogen in side chain
    side_chain_arginine_graph.add_node('C5', atom='C')  
    side_chain_arginine_graph.add_node('N2', atom='N')  
    side_chain_arginine_graph.add_node('N3', atom='N')  
    
    side_chain_arginine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_arginine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_arginine_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond 
    side_chain_arginine_graph.add_edge('C4', 'N', bond_type='1')   # Single bond 
    side_chain_arginine_graph.add_edge('N', 'C5', bond_type='2')   # Double bond 
    side_chain_arginine_graph.add_edge('C5', 'N2', bond_type='1')  # Single bond 
    side_chain_arginine_graph.add_edge('C5', 'N3', bond_type='1')  # Single bond 
    return side_chain_arginine_graph

'''
Creates another graph pattern for arginine side chain (isomer).
'''
def arginine2():
    side_chain_arginine_graph = nx.Graph()
    side_chain_arginine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_arginine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_arginine_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_arginine_graph.add_node('C4', atom='C')  # Carbon in side chain
    side_chain_arginine_graph.add_node('N', atom='N')   # Nitrogen in side chain
    side_chain_arginine_graph.add_node('C5', atom='C')  
    side_chain_arginine_graph.add_node('N2', atom='N')  
    side_chain_arginine_graph.add_node('N3', atom='N') 
    side_chain_arginine_graph.add_node('N4', atom='N')
    
    side_chain_arginine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_arginine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_arginine_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond 
    side_chain_arginine_graph.add_edge('C4', 'N', bond_type='1')   # Single bond 
    side_chain_arginine_graph.add_edge('N', 'C5', bond_type='1')   # Single bond 
    side_chain_arginine_graph.add_edge('C5', 'N2', bond_type='1')  # Single bond 
    side_chain_arginine_graph.add_edge('C5', 'N3', bond_type='2')  # Double bond 
    side_chain_arginine_graph.add_edge('C5', 'N4', bond_type='1')  # Single bond
    return side_chain_arginine_graph

'''
Creates a graph pattern for proline side chain.
'''
def proline():
    side_chain_proline_graph = nx.Graph()
    side_chain_proline_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_proline_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_proline_graph.add_node('C3', atom='C')  # Additional carbon in side chain
    side_chain_proline_graph.add_node('C4', atom='C')  # Additional carbon in side chain
    side_chain_proline_graph.add_node('N', atom='N')   # Nitrogen in side chain (ring closure)
    
    side_chain_proline_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_proline_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_proline_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond
    side_chain_proline_graph.add_edge('C4', 'N', bond_type='1')   # Single bond (connecting to N)
    return side_chain_proline_graph 

'''
Creates a graph pattern for phenylalanine side chain.
'''
def phenylalanine():
    side_chain_phenylalanine_graph = nx.Graph()
    side_chain_phenylalanine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_phenylalanine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_phenylalanine_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_phenylalanine_graph.add_node('C4', atom='C')  # Carbon in side chain
    side_chain_phenylalanine_graph.add_node('C5', atom='C')  # Carbon in side chain
    side_chain_phenylalanine_graph.add_node('C6', atom='C')  # Carbon in side chain (benzene ring)
    side_chain_phenylalanine_graph.add_node('C7', atom='C')  # Carbon in side chain (benzene ring)
    side_chain_phenylalanine_graph.add_node('C8', atom='C')  # Carbon in side chain (benzene ring)
    
    side_chain_phenylalanine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_phenylalanine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_phenylalanine_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond
    side_chain_phenylalanine_graph.add_edge('C3', 'C8', bond_type='2')  # Double bond
    side_chain_phenylalanine_graph.add_edge('C4', 'C5', bond_type='2')  # Double bond
    side_chain_phenylalanine_graph.add_edge('C5', 'C6', bond_type='1')  # Single bond
    side_chain_phenylalanine_graph.add_edge('C6', 'C7', bond_type='2')  # Double bond
    side_chain_phenylalanine_graph.add_edge('C7', 'C8', bond_type='1')  # Single bond

    return side_chain_phenylalanine_graph

'''
Creates another graph pattern for phenylalanine side chain (isomer).
'''
def phenylalanine1():
    side_chain_phenylalanine_graph = nx.Graph()
    side_chain_phenylalanine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_phenylalanine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_phenylalanine_graph.add_node('C3', atom='C')  # Carbon in side chain
    side_chain_phenylalanine_graph.add_node('C4', atom='C')  # Carbon in side chain
    side_chain_phenylalanine_graph.add_node('C5', atom='C')  # Carbon in side chain
    side_chain_phenylalanine_graph.add_node('C6', atom='C')  # Carbon in side chain 
    side_chain_phenylalanine_graph.add_node('C7', atom='C')  # Carbon in side chain 
    side_chain_phenylalanine_graph.add_node('C8', atom='C')  # Carbon in side chain 
    
    side_chain_phenylalanine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_phenylalanine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_phenylalanine_graph.add_edge('C3', 'C4', bond_type='2')  # Double bond
    side_chain_phenylalanine_graph.add_edge('C3', 'C8', bond_type='1')  # Single bond
    side_chain_phenylalanine_graph.add_edge('C4', 'C5', bond_type='1')  # Single bond
    side_chain_phenylalanine_graph.add_edge('C5', 'C6', bond_type='2')  # Double bond
    side_chain_phenylalanine_graph.add_edge('C6', 'C7', bond_type='1')  # Single bond
    side_chain_phenylalanine_graph.add_edge('C7', 'C8', bond_type='2')  # Double bond

    return side_chain_phenylalanine_graph

'''
Creates a graph pattern for tyrosine side chain.
'''
def tyrosine():
    side_chain_tyrosine_graph = nx.Graph()
    side_chain_tyrosine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_tyrosine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_tyrosine_graph.add_node('C3', atom='C')  # Second carbon in side chain
    side_chain_tyrosine_graph.add_node('C4', atom='C')  # Third carbon in side chain
    side_chain_tyrosine_graph.add_node('C5', atom='C')  # Fourth carbon in side chain 
    side_chain_tyrosine_graph.add_node('C6', atom='C')  # Fifth carbon in side chain 
    side_chain_tyrosine_graph.add_node('C7', atom='C')  # Sixth carbon in side chain 
    side_chain_tyrosine_graph.add_node('C8', atom='C')  # Seventh carbon in side chain 
    side_chain_tyrosine_graph.add_node('O', atom='O')   # Hydroxyl group
    
    side_chain_tyrosine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_tyrosine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_tyrosine_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond
    side_chain_tyrosine_graph.add_edge('C3', 'C8', bond_type='2')  # Double bond 
    side_chain_tyrosine_graph.add_edge('C4', 'C5', bond_type='2')  # Double bond 
    side_chain_tyrosine_graph.add_edge('C5', 'C6', bond_type='1')  # Single bond 
    side_chain_tyrosine_graph.add_edge('C6', 'C7', bond_type='2')  # Double bond 
    side_chain_tyrosine_graph.add_edge('C6', 'O', bond_type='1')   # Single bond 
    side_chain_tyrosine_graph.add_edge('C7', 'C8', bond_type='1')  # Single bond 
    return side_chain_tyrosine_graph

'''
Creates another graph pattern for tyrosine side chain (isomer).
'''
def tyrosine1():
    side_chain_tyrosine_graph = nx.Graph()
    side_chain_tyrosine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_tyrosine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_tyrosine_graph.add_node('C3', atom='C')  # Second carbon in side chain
    side_chain_tyrosine_graph.add_node('C4', atom='C')  # Third carbon in side chain
    side_chain_tyrosine_graph.add_node('C5', atom='C')  # Fourth carbon in side chain 
    side_chain_tyrosine_graph.add_node('C6', atom='C')  # Fifth carbon in side chain 
    side_chain_tyrosine_graph.add_node('C7', atom='C')  # Sixth carbon in side chain 
    side_chain_tyrosine_graph.add_node('C8', atom='C')  # Seventh carbon in side chain 
    side_chain_tyrosine_graph.add_node('O', atom='O')   # Hydroxyl group
    
    side_chain_tyrosine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_tyrosine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_tyrosine_graph.add_edge('C3', 'C4', bond_type='2')  # Double bond
    side_chain_tyrosine_graph.add_edge('C3', 'C8', bond_type='1')  # Single bond 
    side_chain_tyrosine_graph.add_edge('C4', 'C5', bond_type='1')  # Single bond 
    side_chain_tyrosine_graph.add_edge('C5', 'C6', bond_type='2')  # Double bond 
    side_chain_tyrosine_graph.add_edge('C6', 'C7', bond_type='1')  # Single bond 
    side_chain_tyrosine_graph.add_edge('C6', 'O', bond_type='1')   # Single bond 
    side_chain_tyrosine_graph.add_edge('C7', 'C8', bond_type='2')  # Double bond 
    return side_chain_tyrosine_graph

'''
Creates a graph pattern for histidine side chain.
'''
def histidine():
    side_chain_histidine_graph = nx.Graph()
    side_chain_histidine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_histidine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_histidine_graph.add_node('C3', atom='C')  # Second carbon in side chain
    side_chain_histidine_graph.add_node('C4', atom='C')  # Third carbon in side chain
    side_chain_histidine_graph.add_node('C5', atom='C')  # Fourth carbon in side chain 
    side_chain_histidine_graph.add_node('N', atom='N')   # Nitrogen in side chain
    side_chain_histidine_graph.add_node('N2', atom='N')  # Second nitrogen in side chain
    
    side_chain_histidine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_histidine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_histidine_graph.add_edge('C3', 'N', bond_type='1')   # Single bond 
    side_chain_histidine_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond
    side_chain_histidine_graph.add_edge('N', 'C5', bond_type='1')   # Single bond 
    side_chain_histidine_graph.add_edge('C5', 'N2', bond_type='2')  # Double bond 
    side_chain_histidine_graph.add_edge('N2', 'C4', bond_type='1')  # Single bond 
    
    return side_chain_histidine_graph

'''
Creates another graph pattern for histidine side chain (isomer).
'''
def histidine1():
    side_chain_histidine_graph = nx.Graph()
    side_chain_histidine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_histidine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_histidine_graph.add_node('C3', atom='C')  # Second carbon in side chain
    side_chain_histidine_graph.add_node('C4', atom='C')  # Third carbon in side chain
    side_chain_histidine_graph.add_node('C5', atom='C')  # Fourth carbon in side chain 
    side_chain_histidine_graph.add_node('N', atom='N')   # Nitrogen in side chain
    side_chain_histidine_graph.add_node('N2', atom='N')  # Second nitrogen in side chain
    
    side_chain_histidine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_histidine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_histidine_graph.add_edge('C3', 'N', bond_type='1')   # Single bond 
    side_chain_histidine_graph.add_edge('C3', 'C4', bond_type='2')  # Double bond
    side_chain_histidine_graph.add_edge('N', 'C5', bond_type='1')   # Single bond
    side_chain_histidine_graph.add_edge('C5', 'N2', bond_type='2')  # Double bond 
    side_chain_histidine_graph.add_edge('N2', 'C4', bond_type='1')  # Single bond 
    
    return side_chain_histidine_graph

'''
Creates another graph pattern for histidine side chain (isomer).
'''
def histidine2():
    side_chain_histidine_graph = nx.Graph()
    side_chain_histidine_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_histidine_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_histidine_graph.add_node('C3', atom='C')  # Second carbon in side chain
    side_chain_histidine_graph.add_node('C4', atom='C')  # Third carbon in side chain
    side_chain_histidine_graph.add_node('C5', atom='C')  # Fourth carbon in side chain 
    side_chain_histidine_graph.add_node('N', atom='N')   # Nitrogen in side chain
    side_chain_histidine_graph.add_node('N2', atom='N')  # Second nitrogen in side chain
    
    side_chain_histidine_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_histidine_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_histidine_graph.add_edge('C3', 'N', bond_type='2')   # Double bond 
    side_chain_histidine_graph.add_edge('C3', 'C4', bond_type='1')  # Single bond
    side_chain_histidine_graph.add_edge('N', 'C5', bond_type='1')   # Single bond 
    side_chain_histidine_graph.add_edge('C5', 'N2', bond_type='2')  # Double bond 
    side_chain_histidine_graph.add_edge('N2', 'C4', bond_type='1')  # Single bond 
    return side_chain_histidine_graph

'''
Creates a graph pattern for tryptophan side chain.
'''
def tryptophan():
    side_chain_tryptophan_graph = nx.Graph()
    side_chain_tryptophan_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_tryptophan_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_tryptophan_graph.add_node('C3', atom='C')  # Second carbon in side chain
    side_chain_tryptophan_graph.add_node('C4', atom='C')  # Third carbon in side chain
    side_chain_tryptophan_graph.add_node('C5', atom='C')  # Fourth carbon in side chain
    side_chain_tryptophan_graph.add_node('C6', atom='C')  # Fifth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C7', atom='C')  # Sixth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C8', atom='C')  # Seventh carbon in side chain
    side_chain_tryptophan_graph.add_node('C9', atom='C')  # Eighth carbon in side chain
    side_chain_tryptophan_graph.add_node('C10', atom='C') # Ninth carbon in side chain 
    side_chain_tryptophan_graph.add_node('N', atom='N')   # Nitrogen in side chain 

    side_chain_tryptophan_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_tryptophan_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_tryptophan_graph.add_edge('C3', 'C4', bond_type='2')  # Double bond
    side_chain_tryptophan_graph.add_edge('C3', 'C5', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C4', 'N', bond_type='1')   # Single bond 
    side_chain_tryptophan_graph.add_edge('N', 'C6', bond_type='1')   # Single bond 
    side_chain_tryptophan_graph.add_edge('C5', 'C6', bond_type='2')  # Double bond 
    side_chain_tryptophan_graph.add_edge('C5', 'C10', bond_type='1') # Single bond 
    side_chain_tryptophan_graph.add_edge('C6', 'C7', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C7', 'C8', bond_type='2')  # Double bond 
    side_chain_tryptophan_graph.add_edge('C8', 'C9', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C9', 'C10', bond_type='2') # Double bond 
    
    return side_chain_tryptophan_graph

'''
Creates another graph pattern for tryptophan side chain (isomer).
'''
def tryptophan1():
    side_chain_tryptophan_graph = nx.Graph()
    side_chain_tryptophan_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_tryptophan_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_tryptophan_graph.add_node('C3', atom='C')  # Second carbon in side chain
    side_chain_tryptophan_graph.add_node('C4', atom='C')  # Third carbon in side chain
    side_chain_tryptophan_graph.add_node('C5', atom='C')  # Fourth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C6', atom='C')  # Fifth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C7', atom='C')  # Sixth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C8', atom='C')  # Seventh carbon in side chain 
    side_chain_tryptophan_graph.add_node('C9', atom='C')  # Eighth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C10', atom='C') # Ninth carbon in side chain 
    side_chain_tryptophan_graph.add_node('N', atom='N')   # Nitrogen in side chain 

    side_chain_tryptophan_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_tryptophan_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_tryptophan_graph.add_edge('C3', 'C4', bond_type='2')  # Double bond
    side_chain_tryptophan_graph.add_edge('C3', 'C5', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C4', 'N', bond_type='1')   # Single bond 
    side_chain_tryptophan_graph.add_edge('N', 'C6', bond_type='2')   # Double bond 
    side_chain_tryptophan_graph.add_edge('C5', 'C6', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C5', 'C10', bond_type='2') # Double bond 
    side_chain_tryptophan_graph.add_edge('C6', 'C7', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C7', 'C8', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C8', 'C9', bond_type='2')  # Double bond 
    side_chain_tryptophan_graph.add_edge('C9', 'C10', bond_type='1') # Single bond 
    
    return side_chain_tryptophan_graph

'''
Creates another graph pattern for tryptophan side chain (isomer).
'''
def tryptophan2():
    side_chain_tryptophan_graph = nx.Graph()
    side_chain_tryptophan_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_tryptophan_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_tryptophan_graph.add_node('C3', atom='C')  # Second carbon in side chain
    side_chain_tryptophan_graph.add_node('C4', atom='C')  # Third carbon in side chain
    side_chain_tryptophan_graph.add_node('C5', atom='C')  # Fourth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C6', atom='C')  # Fifth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C7', atom='C')  # Sixth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C8', atom='C')  # Seventh carbon in side chain
    side_chain_tryptophan_graph.add_node('C9', atom='C')  # Eighth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C10', atom='C') # Ninth carbon in side chain 
    side_chain_tryptophan_graph.add_node('N', atom='N')   # Nitrogen in side chain 

    side_chain_tryptophan_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_tryptophan_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_tryptophan_graph.add_edge('C3', 'C4', bond_type='2')  # Double bond
    side_chain_tryptophan_graph.add_edge('C3', 'C5', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C4', 'N', bond_type='1')   # Single bond 
    side_chain_tryptophan_graph.add_edge('N', 'C6', bond_type='1')   # Single bond 
    side_chain_tryptophan_graph.add_edge('C5', 'C6', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C5', 'C10', bond_type='2') # Double bond 
    side_chain_tryptophan_graph.add_edge('C6', 'C7', bond_type='2')  # Double bond 
    side_chain_tryptophan_graph.add_edge('C7', 'C8', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C8', 'C9', bond_type='2')  # Double bond 
    side_chain_tryptophan_graph.add_edge('C9', 'C10', bond_type='1') # Single bond 
    
    return side_chain_tryptophan_graph

'''
Creates another graph pattern for tryptophan side chain (isomer).
'''
def tryptophan3():
    side_chain_tryptophan_graph = nx.Graph()
    side_chain_tryptophan_graph.add_node('Cα', atom='C')  # Central carbon (Cα)
    side_chain_tryptophan_graph.add_node('C2', atom='C')  # Carbon in side chain
    side_chain_tryptophan_graph.add_node('C3', atom='C')  # Second carbon in side chain
    side_chain_tryptophan_graph.add_node('C4', atom='C')  # Third carbon in side chain
    side_chain_tryptophan_graph.add_node('C5', atom='C')  # Fourth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C6', atom='C')  # Fifth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C7', atom='C')  # Sixth carbon in side chain 
    side_chain_tryptophan_graph.add_node('C8', atom='C')  # Seventh carbon in side chain 
    side_chain_tryptophan_graph.add_node('C9', atom='C')  # Eighth carbon in side chain
    side_chain_tryptophan_graph.add_node('C10', atom='C') # Ninth carbon in side chain 
    side_chain_tryptophan_graph.add_node('N', atom='N')   # Nitrogen in side chain 

    side_chain_tryptophan_graph.add_edge('Cα', 'C2', bond_type='1')  # Single bond
    side_chain_tryptophan_graph.add_edge('C2', 'C3', bond_type='1')  # Single bond
    side_chain_tryptophan_graph.add_edge('C3', 'C4', bond_type='2')  # Double bond
    side_chain_tryptophan_graph.add_edge('C3', 'C5', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C4', 'N', bond_type='1')   # Single bond 
    side_chain_tryptophan_graph.add_edge('N', 'C6', bond_type='2')   # Double bond 
    side_chain_tryptophan_graph.add_edge('C5', 'C6', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C5', 'C10', bond_type='1') # Single bond 
    side_chain_tryptophan_graph.add_edge('C6', 'C7', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C7', 'C8', bond_type='2')  # Double bond 
    side_chain_tryptophan_graph.add_edge('C8', 'C9', bond_type='1')  # Single bond 
    side_chain_tryptophan_graph.add_edge('C9', 'C10', bond_type='2') # Double bond 
    return side_chain_tryptophan_graph

# Dictionary of all amino acid patterns
amino_acids = {'ALA': alanine(), 'GLY': glycine(),'VAL': valine(),'LEU': leucine(),'ILE': isoleucine(),'PRO': proline(),'SER': serine(),
               'MET': methionine(),'CYS': cysteine(), 'PHE': phenylalanine(), 'PHE1': phenylalanine1(),'ASN': asparagine(),
               'TYR': tyrosine(), 'TYR1':tyrosine1(), 'ASP': aspartate(), 'HIS': histidine(),'HIS1': histidine1(),'HIS2': histidine2(),
               'TRP': tryptophan(),'TRP1':tryptophan1(), 'TRP2': tryptophan2(),'TRP3': tryptophan3(),
               'ARG': arginine(), 'ARG1': arginine1(), 'ARG2': arginine2(),
               'GLU': glutamate(), 'THR': threonine(), 'GLN': glutamine(), 'LYS': lysine()}

'''
Compares two graph structures for isomorphism.
First checks if node and edge counts match, then compares atom attributes.
Finally, checks if C-alpha atoms correspond between graphs.
Returns True if graphs are isomorphic, False otherwise.
'''
def compare_graphs(graph1, graph2):
    if len(graph1.nodes) != len(graph2.nodes) or len(graph1.edges) != len(graph2.edges):
        return False
    # Compare nodes by their 'atom' attributes
    nodes1 = sorted(attr for _, attr in graph1.nodes(data='atom'))
    nodes2 = sorted(attr for _, attr in graph2.nodes(data='atom'))
    if nodes1 != nodes2:
        return False
    root1 = find_c_alpha(graph1, find_all_backbones_and_C_alpha(Graf)[1])
    root2 = 'Cα' 
    node_match = lambda x, y: x['atom'] == y['atom']
    edge_match = lambda x, y: x['bond_type'] == y['bond_type']

    GM = GraphMatcher(graph1, graph2, node_match=node_match, edge_match=edge_match)

    # Check isomorphisms where root1 maps to 'Cα'
    for mapping in GM.isomorphisms_iter():
        if mapping.get(root1) == root2:
            return True
    return False

'''
Finds the C-alpha node in a graph using a list of known C-alpha nodes.
Returns the C-alpha node if found, None otherwise.
'''
def find_c_alpha(graph, c_alphas):
    for node in graph.nodes():
        if graph.nodes[node]['atom'] == 'C' and node in c_alphas:
            return node
    return None

'''
Identifies an amino acid by comparing its side chain graph to known patterns.
Returns the amino acid code if identified, None otherwise.
'''
def identify_ac(graph, amino_acids, c_alphas):
    result = None    
    for name in amino_acids.keys():
        side_chain = amino_acids[name]
        c_alpha = find_c_alpha(graph, c_alphas)
        if compare_graphs(graph, side_chain):
            result = name
            for letter in result:
                if letter.isnumeric():
                    result = result.replace(letter, '')
            
            break
    return result

'''
Converts a NetworkX graph to the format required by GraKeL for graph kernel computations.
Returns a tuple containing edges, node labels, and edge labels.
'''
def nx_to_grakel(nx_graph):
    # Map node names to integers (GraKeL expects this)
    node_mapping = {n: i for i, n in enumerate(nx_graph.nodes())}
    # Node labels
    node_labels = {node_mapping[n]: nx_graph.nodes[n].get('atom', '') for n in nx_graph.nodes()}
    # Edges with labels (GraKeL expects edges as (u, v), not (u, v, attr))
    edges = [(node_mapping[u], node_mapping[v]) for u, v in nx_graph.edges()]
    edge_labels = {(node_mapping[u], node_mapping[v]): nx_graph[u][v].get('bond_type', '')
                   for u, v in nx_graph.edges()}
    return (edges, node_labels, edge_labels)

'''
Calculates the similarity between two graphs using the Weisfeiler-Lehman graph kernel.
Returns a similarity score between 0 and 1, where 1 means perfect match.
'''
def wl_similarity(graph1, graph2):
    # Convert to GraKel format
    gk_graph1 = nx_to_grakel(graph1)
    gk_graph2 = nx_to_grakel(graph2)
    
    # Initialize and compute the kernel
    kernel = GraphKernel(
        kernel=[{"name": "weisfeiler_lehman", "n_iter": 3},{"name": "subtree_wl"}],  # 3 iterations
        normalize=True  # Normalize scores to [0, 1]
    )
    K = kernel.fit_transform([gk_graph1, gk_graph2])
    return K[0][1]

'''
Finds the amino acid pattern with the highest similarity to a given side chain.
Returns a tuple with the amino acid name and its similarity score.
'''
def similarity(side_ch, backbone):
    score = 0
    res = None
    for name in backbone.keys():
        ref = backbone[name]
        try:
            if wl_similarity(side_ch, ref) > score:
                score = wl_similarity(side_ch, ref)
                res = name
                
        except Exception as e:
            continue
    return res, score


'''
Maps each atom to its corresponding amino acid based on the side chain analysis.
Returns a dictionary where keys are atom IDs and values are amino acid codes.
'''
def map_atoms_to_amino_acids(data):
    Graf = get_nodes_and_edges(data)
    backbone = find_all_backbones_and_C_alpha(Graf)[0]
    c_alphas = find_all_backbones_and_C_alpha(Graf)[1]
    atom_to_aa = {}

    # First identify all side chains and their amino acid types
    side_chain_atoms = {}
    for side in c_alphas:
        side_chain = extract_side_chain(side, backbone, c_alphas, Graf)
        ac = identify_ac(side_chain, amino_acids, c_alphas)

        if not ac:
            if similarity(side_chain, amino_acids)[1] > 0.85:
                ac = similarity(side_chain, amino_acids)[0]
            else:
                ac = "UNK"  # Unknown amino acid
        else:
            # Remove numeric identifiers from amino acid names
            ac = ''.join([letter for letter in ac if not letter.isnumeric()])

        # Map all atoms in this side chain to this amino acid
        for atom in side_chain.nodes():
            side_chain_atoms[atom] = ac

    # Now map backbone atoms to their corresponding amino acids
    # We'll assign backbone atoms to the nearest C-alpha
    c_alpha_to_aa = {}
    for side in c_alphas:
        if side in side_chain_atoms:
            c_alpha_to_aa[side] = side_chain_atoms[side]

    # For each backbone atom, find the closest C-alpha
    for atom in backbone:
        if atom in side_chain_atoms:
            continue  # Already assigned

        closest_ca = None
        min_distance = float('inf')

        # Simple distance calculation (in graph space)
        for ca in c_alphas:
            try:
                distance = nx.shortest_path_length(Graf, atom, ca)
                if distance < min_distance:
                    min_distance = distance
                    closest_ca = ca
            except nx.NetworkXNoPath:
                continue

        if closest_ca and closest_ca in c_alpha_to_aa:
            atom_to_aa[atom] = c_alpha_to_aa[closest_ca]
        else:
            atom_to_aa[atom] = "UNK"  # Default if we can't assign

    # Combine the mappings
    atom_to_aa.update(side_chain_atoms)

    return atom_to_aa

'''
Creates a new SDF file with amino acid annotations for each atom.
The new file contains the original data plus amino acid labels.
'''
def create_annotated_sdf(input_data, output_file):
    # Get atom to amino acid mapping
    atom_to_aa = map_atoms_to_amino_acids(input_data)

    # Write the new SDF file
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
                    atom_id = parts[2]
                    if atom_id in atom_to_aa:
                        # Keep original line and append amino acid
                        outfile.write(f"{line.rstrip()} {atom_to_aa[atom_id]}\n")
                    else:
                        outfile.write(f"{line.rstrip()} UNK\n")
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

# Main execution
Graf = get_nodes_and_edges(raw_data)

# Create the annotated SDF file
import os
base_name = os.path.splitext(os.path.basename(file_path))[0]
output_sdf_path = f"{base_name}_annotated.sdf"
print(f"\n\nCreating annotated SDF file at {output_sdf_path}")
create_annotated_sdf(raw_data, output_sdf_path)
print(f"Finished creating annotated SDF file with amino acid annotations.")
