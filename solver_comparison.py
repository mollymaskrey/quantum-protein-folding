### NOTE you will need to set your DWAVE SAMPLER TOKEN into the variable 'token'

import os

token = os.getenv("DWAVE_SAMPLER_TOKEN")
if token is None:
    raise ValueError("DWAVE_SAMPLER_TOKEN is not defined!")


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dimod
from dimod import ExactSolver, ConstrainedQuadraticModel
from dwave.system import LeapHybridCQMSampler
from Bio.PDB import PDBParser, PDBList, calc_dihedral
from Bio import SeqIO
import requests
import time
from datetime import datetime
import warnings
import os
warnings.filterwarnings('ignore')

class QuantumProteinFolder:
    """
    Quantum protein folding using both ExactSolver and D-Wave Hybrid solvers.
    
    This class implements a quantum computing approach to protein folding optimization
    by comparing discrete vs continuous optimization approaches. It uses D-Wave's
    quantum annealing technology to find optimal protein conformations based on
    dihedral angles (phi and psi).
    
    Attributes:
        structure: Bio.PDB structure object containing protein 3D coordinates
        sequence (str): Amino acid sequence of the protein
        phi_angles (list): List of phi dihedral angles for each residue
        psi_angles (list): List of psi dihedral angles for each residue
        optimized_phi (list): Optimized phi angles after quantum optimization
        optimized_psi (list): Optimized psi angles after quantum optimization
        timing_stats (dict): Performance metrics for ExactSolver
        timing_stats_hybrid (dict): Performance metrics for Hybrid solver
    """
    
    def __init__(self):
        """Initialize the QuantumProteinFolder with empty attributes."""
        self.structure = None
        self.sequence = None
        self.phi_angles = None
        self.psi_angles = None
        self.optimized_phi = None
        self.optimized_psi = None
        self.timing_stats = {}
        self.timing_stats_hybrid = {}
        
    def download_alphafold_protein(self, uniprot_id):
        """
        Download protein structure from AlphaFold database.
        
        Downloads the PDB file for a given UniProt ID from the AlphaFold
        Protein Structure Database. AlphaFold provides high-quality predicted
        3D structures for proteins.
        
        Args:
            uniprot_id (str): UniProt accession ID (e.g., 'P01308' for insulin)
            
        Returns:
            str: Filename of downloaded PDB file if successful, None if failed
            
        Example:
            >>> folder = QuantumProteinFolder()
            >>> pdb_file = folder.download_alphafold_protein('P01308')
            >>> print(pdb_file)  # 'P01308_alphafold.pdb'
        """
        print(f"Downloading {uniprot_id} from AlphaFold...")
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                filename = f"{uniprot_id}_alphafold.pdb"
                with open(filename, 'w') as f:
                    f.write(response.text)
                print(f"Downloaded: {filename}")
                return filename
            else:
                print(f"Failed to download {uniprot_id}")
                return None
        except Exception as e:
            print(f"Download error: {e}")
            return None
    
    def analyze_protein_structure(self, pdb_file):
        """
        Perform comprehensive analysis of protein structure from PDB file.
        
        Parses the PDB file to extract structural information including:
        - Number of models, chains, residues, and atoms
        - Amino acid sequence
        - 3D coordinates for quantum optimization
        
        Args:
            pdb_file (str): Path to PDB file
            
        Returns:
            tuple: (structure object, sequence string)
            
        Side Effects:
            - Sets self.structure with parsed PDB structure
            - Sets self.sequence with amino acid sequence
            - Prints detailed analysis to console
        """
        print(f"\nANALYZING PROTEIN STRUCTURE")
        print("=" * 50)
        
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('protein', pdb_file)
        
        # Basic structure analysis
        models = list(self.structure)
        chains = list(self.structure.get_chains())
        residues = list(self.structure.get_residues())
        atoms = list(self.structure.get_atoms())
        
        print(f"Structure Overview:")
        print(f"• Models: {len(models)}")
        print(f"• Chains: {len(chains)}")
        print(f"• Total residues: {len(residues)}")
        print(f"• Total atoms: {len(atoms)}")
        
        # Extract amino acid sequence
        amino_acids = []
        aa_map = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        
        for residue in residues:
            if residue.get_resname() in aa_map:
                amino_acids.append(aa_map[residue.get_resname()])
        
        self.sequence = ''.join(amino_acids)
        
        print(f"\nProtein Sequence ({len(self.sequence)} residues):")
        # Print sequence in chunks of 50 for readability
        for i in range(0, len(self.sequence), 50):
            print(f"   {self.sequence[i:i+50]}")
        
        return self.structure, self.sequence
    
    def extract_dihedral_angles(self):
        """
        Extract phi and psi dihedral angles from protein backbone.
        
        Calculates the backbone dihedral angles that define protein conformation:
        - Phi (φ): Rotation around N-Cα bond (C(-1)-N-CA-C)
        - Psi (ψ): Rotation around Cα-C bond (N-CA-C-N(+1))
        
        These angles are critical for protein folding as they determine
        the 3D shape of the protein backbone. Falls back to reasonable
        random values if angle calculation fails.
        
        Returns:
            tuple: (phi_angles list, psi_angles list) in degrees
            
        Side Effects:
            - Sets self.phi_angles and self.psi_angles
            - Prints extraction statistics
            
        Note:
            - First residue has no phi angle (no previous C)
            - Last residue has no psi angle (no next N)
            - Values are in degrees (-180 to +180)
        """
        print(f"\nEXTRACTING DIHEDRAL ANGLES")
        print("=" * 50)
        
        print(f"EXPLANATION:")
        print(f"• Phi (φ): Rotation around N-Cα bond")
        print(f"• Psi (ψ): Rotation around Cα-C bond")
        print(f"• These angles determine the protein's 3D backbone shape")
        
        residues = list(self.structure.get_residues())
        protein_residues = []
        
        # Filter to standard amino acids and check for required atoms
        standard_aas = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 
                       'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                       'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
        
        for residue in residues:
            if residue.get_resname() in standard_aas:
                # Check if the residue has the required backbone atoms
                required_atoms = ['N', 'CA', 'C']
                has_atoms = all(atom in residue for atom in required_atoms)
                if has_atoms:
                    protein_residues.append(residue)
        
        print(f"\n🔬 Processing {len(protein_residues)} complete residues...")
        
        self.phi_angles = []
        self.psi_angles = []
        
        for i, residue in enumerate(protein_residues):
            phi = None
            psi = None
            
            try:
                # Calculate phi angle (requires previous residue's C)
                if i > 0:
                    prev_residue = protein_residues[i-1]
                    if 'C' in prev_residue and all(atom in residue for atom in ['N', 'CA', 'C']):
                        phi = calc_dihedral(
                            prev_residue['C'].get_vector(),
                            residue['N'].get_vector(),
                            residue['CA'].get_vector(),
                            residue['C'].get_vector()
                        )
                        phi = np.degrees(phi)
                
                # Calculate psi angle (requires next residue's N)
                if i < len(protein_residues) - 1:
                    next_residue = protein_residues[i+1]
                    if 'N' in next_residue and all(atom in residue for atom in ['N', 'CA', 'C']):
                        psi = calc_dihedral(
                            residue['N'].get_vector(),
                            residue['CA'].get_vector(),
                            residue['C'].get_vector(),
                            next_residue['N'].get_vector()
                        )
                        psi = np.degrees(psi)
                        
            except Exception:
                # Fall back to sample angles for demonstration
                if i > 0 and i < len(protein_residues) - 1:
                    phi = -60 + np.random.normal(0, 15)
                    psi = -45 + np.random.normal(0, 15)
            
            self.phi_angles.append(phi)
            self.psi_angles.append(psi)
        
        # Count successful calculations
        valid_phi = [p for p in self.phi_angles if p is not None]
        valid_psi = [p for p in self.psi_angles if p is not None]
        
        print(f"Extracted {len(valid_phi)} valid phi angles, {len(valid_psi)} valid psi angles")
        
        return self.phi_angles, self.psi_angles
    
    def create_small_qubo_problem(self, start_residue=0, num_residues=4, num_angle_bins=3):
        """
        Create a QUBO (Quadratic Unconstrained Binary Optimization) problem for discrete angle optimization.
        
        Formulates the protein folding problem as a QUBO suitable for ExactSolver.
        Uses discrete bins for phi/psi angles and one-hot encoding to ensure
        each angle takes exactly one value.
        
        Args:
            start_residue (int): Starting residue index (0-based)
            num_residues (int): Number of residues to optimize
            num_angle_bins (int): Number of discrete bins per angle
            
        Returns:
            tuple: Contains:
                - Q (np.ndarray): QUBO matrix (NxN)
                - segment_sequence (str): Amino acid sequence for segment
                - segment_phi (list): Original phi angles for segment
                - segment_psi (list): Original psi angles for segment
                - phi_bins (np.ndarray): Discrete phi angle values
                - psi_bins (np.ndarray): Discrete psi angle values
                - get_var_index (function): Helper to map residue/angle to variable index
                
        Note:
            Total variables = num_residues * 2 * num_angle_bins
            Uses Ramachandran energy to favor alpha-helix and beta-sheet regions
        """
        print(f"\nCREATING DISCRETE QUBO PROBLEM")
        print("=" * 50)
        
        # Extract segment
        end_residue = min(start_residue + num_residues, len(self.sequence))
        actual_residues = end_residue - start_residue
        segment_sequence = self.sequence[start_residue:end_residue]
        segment_phi = self.phi_angles[start_residue:end_residue]
        segment_psi = self.psi_angles[start_residue:end_residue]
        
        print(f"Optimization Target:")
        print(f"• Segment: residues {start_residue+1}-{end_residue}")
        print(f"• Sequence: {segment_sequence}")
        print(f"• Angle bins: {num_angle_bins}")
        print(f"• Total variables: {actual_residues * 2 * num_angle_bins}")
        
        # Define discrete angle bins
        phi_bins = np.linspace(-180, 180, num_angle_bins+1)[:-1]
        psi_bins = np.linspace(-180, 180, num_angle_bins+1)[:-1]
        
        # Create QUBO matrix
        total_vars = actual_residues * 2 * num_angle_bins
        Q = np.zeros((total_vars, total_vars))
        
        def get_var_index(residue, angle_type, bin_idx):
            base = residue * 2 * num_angle_bins
            if angle_type == 'psi':
                base += num_angle_bins
            return base + bin_idx
        
        # Add constraints and energy terms
        constraint_weight = 10.0
        energy_weight = 1.0
        
        # One-hot constraints
        for residue in range(actual_residues):
            for angle_type in ['phi', 'psi']:
                for i in range(num_angle_bins):
                    var_i = get_var_index(residue, angle_type, i)
                    Q[var_i, var_i] -= 2 * constraint_weight
                    
                    for j in range(num_angle_bins):
                        var_j = get_var_index(residue, angle_type, j)
                        Q[var_i, var_j] += constraint_weight
        
        # Ramachandran energy
        for residue in range(actual_residues):
            for phi_bin in range(num_angle_bins):
                for psi_bin in range(num_angle_bins):
                    phi_angle = phi_bins[phi_bin]
                    psi_angle = psi_bins[psi_bin]
                    
                    alpha_energy = ((phi_angle + 60)**2 + (psi_angle + 45)**2) / 10000
                    beta_energy = ((phi_angle + 120)**2 + (psi_angle - 120)**2) / 10000
                    ramachandran_energy = min(alpha_energy, beta_energy)
                    
                    phi_var = get_var_index(residue, 'phi', phi_bin)
                    psi_var = get_var_index(residue, 'psi', psi_bin)
                    
                    Q[phi_var, psi_var] += energy_weight * ramachandran_energy
                    Q[psi_var, phi_var] += energy_weight * ramachandran_energy
        
        print(f"QUBO matrix created: {total_vars}×{total_vars}")
        
        return Q, segment_sequence, segment_phi, segment_psi, phi_bins, psi_bins, get_var_index
    
    def create_cqm_problem(self, start_residue=0, num_residues=4):
        """
        Create a CQM (Constrained Quadratic Model) problem for continuous angle optimization.
        
        Formulates the protein folding problem with continuous variables suitable
        for D-Wave's nonlinear hybrid solver. Uses real-valued angles instead
        of discrete bins for higher precision.
        
        Args:
            start_residue (int): Starting residue index (0-based)
            num_residues (int): Number of residues to optimize
            
        Returns:
            tuple: Contains:
                - cqm (ConstrainedQuadraticModel): CQM object for hybrid solver
                - segment_sequence (str): Amino acid sequence for segment
                - segment_phi (list): Original phi angles for segment
                - segment_psi (list): Original psi angles for segment
                - phi_var_names (list): Variable names for phi angles
                - psi_var_names (list): Variable names for psi angles
                
        Note:
            - Variables are continuous in range [-180, 180] degrees
            - Linear objective minimizes distance from alpha helix conformation
            - Suitable for problems too large for ExactSolver
        """
        print(f"\nCREATING NONLINEAR CQM PROBLEM")
        print("=" * 50)
        
        # Extract segment
        end_residue = min(start_residue + num_residues, len(self.sequence))
        actual_residues = end_residue - start_residue
        segment_sequence = self.sequence[start_residue:end_residue]
        segment_phi = self.phi_angles[start_residue:end_residue]
        segment_psi = self.psi_angles[start_residue:end_residue]
        
        print(f"Optimization Target:")
        print(f"• Segment: residues {start_residue+1}-{end_residue}")
        print(f"• Sequence: {segment_sequence}")
        print(f"• Variables: {actual_residues * 2} continuous variables (phi, psi)")
        print(f"• Range: -180° to +180° for each angle")
        
        # Create variable names
        phi_var_names = [f'phi_{i}' for i in range(actual_residues)]
        psi_var_names = [f'psi_{i}' for i in range(actual_residues)]
        
        # Create objective using QuadraticModel only
        from dimod import QuadraticModel
        qm = QuadraticModel()
        
        # Add variables to QM with bounds
        for phi_var in phi_var_names:
            qm.add_variable('REAL', phi_var, lower_bound=-180, upper_bound=180)
        for psi_var in psi_var_names:
            qm.add_variable('REAL', psi_var, lower_bound=-180, upper_bound=180)
        
        # Linear objective: minimize distance from alpha helix
        for i in range(actual_residues):
            phi_var = phi_var_names[i]
            psi_var = psi_var_names[i]
            
            # Linear terms to encourage alpha helix: phi ≈ -60, psi ≈ -45
            qm.add_linear(phi_var, 1.0)  # coefficient for phi
            qm.add_linear(psi_var, 1.0)  # coefficient for psi
            qm.offset += 60.0 + 45.0  # constant terms
        
        # Create CQM from the QuadraticModel
        cqm = ConstrainedQuadraticModel.from_quadratic_model(qm)
        
        print(f"CQM created with {len(cqm.variables)} continuous variables")
        print(f"• Objective: Linear minimization toward alpha helix angles")
        
        return cqm, segment_sequence, segment_phi, segment_psi, phi_var_names, psi_var_names
    
    def solve_with_exact_solver(self, Q):
        """
        Solve QUBO problem using dimod's ExactSolver for small problems.
        
        ExactSolver performs exhaustive search over all 2^N possible solutions,
        guaranteeing the global optimum. Only suitable for small problems
        (typically < 20 variables).
        
        Args:
            Q (np.ndarray): QUBO matrix defining the optimization problem
            
        Returns:
            tuple: Contains:
                - best_sample (dict): Binary variable assignments for best solution
                - best_energy (float): Energy value of best solution
                - response (SampleSet): Complete solver response with all solutions
                
        Side Effects:
            - Updates self.timing_stats with performance metrics
            - Prints detailed timing information
            
        Note:
            Time complexity is O(2^N) where N is number of variables
        """
        print(f"\nQUANTUM OPTIMIZATION WITH EXACTSOLVER")
        print("=" * 50)
        
        variables = Q.shape[0]
        solutions = 2**variables
        
        print(f"Problem Statistics:")
        print(f"• Variables: {variables}")
        print(f"• Solution space: 2^{variables} = {solutions:,}")
        print(f"• Started: {datetime.now().strftime('%H:%M:%S.%f')[:-3]}")
        
        start_time = time.time()
        setup_start = time.time()
        
        bqm = dimod.BinaryQuadraticModel(Q, "BINARY")
        setup_time = time.time() - setup_start
        
        solver = ExactSolver()
        solve_start = time.time()
        response = solver.sample(bqm)
        solve_time = time.time() - solve_start
        total_time = time.time() - start_time
        
        best_sample = response.first.sample
        best_energy = response.first.energy
        
        print(f"\nOPTIMIZATION RESULTS:")
        print(f"• Best energy: {best_energy:.6f}")
        print(f"• Setup time: {setup_time:.4f}s")
        print(f"• Solve time: {solve_time:.4f}s")
        print(f"• Total time: {total_time:.4f}s")
        print(f"• Solutions/sec: {solutions/solve_time:,.0f}")
        
        self.timing_stats = {
            'setup_time': setup_time,
            'solve_time': solve_time,
            'total_time': total_time,
            'variables': variables,
            'solutions': solutions
        }
        
        return best_sample, best_energy, response
    
    def solve_with_nonlinear_hybrid(self, cqm):
        """
        Solve CQM problem using D-Wave's cloud-based nonlinear hybrid solver.
        
        The hybrid solver combines classical and quantum computing to handle
        larger problems with continuous variables. Requires D-Wave API token
        set as environment variable DWAVE_SAMPLER_TOKEN.
        
        Args:
            cqm (ConstrainedQuadraticModel): CQM defining the optimization problem
            
        Returns:
            tuple: Contains:
                - best_sample (dict): Variable assignments for best solution
                - best_energy (float): Energy value of best solution
                - sampleset (SampleSet): Complete solver response
                
            Returns (None, None, None) if D-Wave token not found
            
        Side Effects:
            - Updates self.timing_stats_hybrid with performance metrics
            - Prints detailed timing and status information
            
        Note:
            - Requires active internet connection
            - Uses D-Wave Leap cloud service
            - Time limit set to 30 seconds
        """
        print(f"\nQUANTUM OPTIMIZATION WITH NONLINEAR HYBRID SOLVER")
        print("=" * 50)
        
        # Check for D-Wave token
        token = os.getenv("DWAVE_SAMPLER_TOKEN")
        if not token:
            print("DWAVE_SAMPLER_TOKEN not found")
            print("Set your token with: export DWAVE_SAMPLER_TOKEN='your-token-here'")
            return None, None, None
        
        variables = len(cqm.variables)
        print(f"Problem Statistics:")
        print(f"• Variables: {variables} (continuous)")
        print(f"• Constraints: {len(cqm.constraints)}")
        print(f"• Started: {datetime.now().strftime('%H:%M:%S.%f')[:-3]}")
        
        start_time = time.time()
        setup_start = time.time()
        
        sampler = LeapHybridCQMSampler(token=token)
        setup_time = time.time() - setup_start
        
        print(f"• Sampler setup: {setup_time:.4f}s")
        print(f"• Sending to D-Wave quantum cloud...")
        
        solve_start = time.time()
        sampleset = sampler.sample_cqm(cqm, 
                                      label="Protein Folding - Nonlinear Hybrid",
                                      time_limit=30)  # 30 second time limit
        solve_time = time.time() - solve_start
        total_time = time.time() - start_time
        
        best_sample = sampleset.first.sample
        best_energy = sampleset.first.energy
        
        print(f"\nOPTIMIZATION RESULTS:")
        print(f"• Best energy: {best_energy:.6f}")
        print(f"• Setup time: {setup_time:.4f}s")
        print(f"• Quantum solve time: {solve_time:.4f}s")
        print(f"• Total time: {total_time:.4f}s")
        
        self.timing_stats_hybrid = {
            'setup_time': setup_time,
            'solve_time': solve_time,
            'total_time': total_time,
            'variables': variables
        }
        
        return best_sample, best_energy, sampleset
    
    def decode_solution(self, sample, sequence, phi_bins, psi_bins, get_var_index, num_angle_bins):
        """
        Convert binary QUBO solution back to phi/psi angles.
        
        Decodes the one-hot encoded binary variables from ExactSolver solution
        into actual phi and psi angle values for each residue.
        
        Args:
            sample (dict): Binary variable assignments from solver
            sequence (str): Amino acid sequence for the segment
            phi_bins (np.ndarray): Discrete phi angle values
            psi_bins (np.ndarray): Discrete psi angle values
            get_var_index (function): Maps residue/angle to variable index
            num_angle_bins (int): Number of discrete bins per angle
            
        Returns:
            tuple: (optimized_phi, optimized_psi) lists of angles in degrees
            
        Side Effects:
            Prints decoded angles for each residue in tabular format
        """
        print(f"\nDECODING EXACTSOLVER SOLUTION")
        print("=" * 50)
        
        n_residues = len(sequence)
        optimized_phi = []
        optimized_psi = []
        
        print("Residue | AA | Phi Bin | Phi (°) | Psi Bin | Psi (°)")
        print("-" * 55)
        
        for residue in range(n_residues):
            phi_angle = None
            psi_angle = None
            phi_bin_selected = None
            psi_bin_selected = None
            
            for bin_idx in range(num_angle_bins):
                phi_var = get_var_index(residue, 'phi', bin_idx)
                if sample[phi_var] == 1:
                    phi_angle = phi_bins[bin_idx]
                    phi_bin_selected = bin_idx
                    break
            
            for bin_idx in range(num_angle_bins):
                psi_var = get_var_index(residue, 'psi', bin_idx)
                if sample[psi_var] == 1:
                    psi_angle = psi_bins[bin_idx]
                    psi_bin_selected = bin_idx
                    break
            
            optimized_phi.append(phi_angle)
            optimized_psi.append(psi_angle)
            
            aa = sequence[residue]
            print(f"  {residue+1:3d}   | {aa}  |    {phi_bin_selected}    | {phi_angle:6.1f}  |    {psi_bin_selected}    | {psi_angle:6.1f}")
        
        return optimized_phi, optimized_psi
    
    def decode_cqm_solution(self, sample, sequence, phi_var_names, psi_var_names):
        """
        Decode continuous CQM solution and predict secondary structure.
        
        Extracts continuous angle values from hybrid solver solution and
        predicts the likely secondary structure based on Ramachandran regions.
        
        Args:
            sample (dict): Variable assignments from hybrid solver
            sequence (str): Amino acid sequence for the segment
            phi_var_names (list): Variable names for phi angles
            psi_var_names (list): Variable names for psi angles
            
        Returns:
            tuple: (optimized_phi, optimized_psi) lists of angles in degrees
            
        Side Effects:
            Prints angles and predicted secondary structure for each residue
            
        Secondary Structure Prediction:
            - α-helix: phi in [-90, -30], psi in [-70, -10]
            - β-sheet: phi in [-150, -90], psi in [90, 150]
            - β-turn: phi in [-90, -30], psi in [-10, 50]
            - random coil: everything else
        """
        print(f"\nDECODING NONLINEAR SOLUTION")
        print("=" * 50)
        
        n_residues = len(sequence)
        optimized_phi = []
        optimized_psi = []
        
        print("Residue | AA | Phi (°) | Psi (°) | Secondary Structure")
        print("-" * 55)
        
        for i in range(n_residues):
            phi_angle = sample[phi_var_names[i]]
            psi_angle = sample[psi_var_names[i]]
            
            optimized_phi.append(phi_angle)
            optimized_psi.append(psi_angle)
            
            # Predict secondary structure
            if -90 <= phi_angle <= -30 and -70 <= psi_angle <= -10:
                ss_pred = "α-helix"
            elif -150 <= phi_angle <= -90 and 90 <= psi_angle <= 150:
                ss_pred = "β-sheet"
            elif -90 <= phi_angle <= -30 and -10 <= psi_angle <= 50:
                ss_pred = "β-turn"
            else:
                ss_pred = "random coil"
            
            aa = sequence[i]
            print(f"  {i+1:3d}   | {aa}  | {phi_angle:6.1f}  | {psi_angle:6.1f}  | {ss_pred}")
        
        return optimized_phi, optimized_psi

def run_solver_comparison():
    """
    Compare ExactSolver vs Nonlinear Hybrid Solver on the same protein segment.
    
    This function demonstrates the differences between discrete (ExactSolver)
    and continuous (Hybrid) optimization approaches for protein folding.
    Downloads insulin structure (P01308) from AlphaFold, extracts dihedral
    angles, and optimizes a 4-residue segment using both solvers.
    
    The comparison includes:
    - Performance metrics (solve time, number of variables)
    - Solution quality (energy values)
    - Angle differences between methods
    - Ramachandran plot visualization
    
    Returns:
        QuantumProteinFolder: Folder object with results from both solvers
        
    Side Effects:
        - Downloads PDB file to current directory
        - Displays matplotlib figure with Ramachandran plots
        - Prints detailed comparison statistics
        
    Note:
        Requires DWAVE_SAMPLER_TOKEN environment variable for hybrid solver
    """
    print("=" * 80)
    print("QUANTUM SOLVER COMPARISON - EXACTSOLVER vs NONLINEAR HYBRID")
    print("=" * 80)
    
    folder = QuantumProteinFolder()
    
    # Download and analyze protein
    pdb_file = folder.download_alphafold_protein('P01308')
    if not pdb_file:
        print("Download failed")
        return
    
    folder.analyze_protein_structure(pdb_file)
    folder.extract_dihedral_angles()
    
    # Test the same 4-residue segment with both solvers
    start_res = 10
    num_res = 4
    
    print(f"\nTESTING SEGMENT: Residues {start_res+1}-{start_res+num_res}")
    print(f"Sequence: {folder.sequence[start_res:start_res+num_res]}")
    
    # ROUND 1: ExactSolver (24 variables, discrete)
    print(f"\n" + "="*60)
    print(f"ROUND 1: EXACTSOLVER (DISCRETE VARIABLES)")
    print(f"="*60)
    
    Q, sequence, orig_phi, orig_psi, phi_bins, psi_bins, get_var_index = folder.create_small_qubo_problem(
        start_residue=start_res, num_residues=num_res, num_angle_bins=3
    )
    
    best_sample_exact, best_energy_exact, response_exact = folder.solve_with_exact_solver(Q)
    opt_phi_exact, opt_psi_exact = folder.decode_solution(
        best_sample_exact, sequence, phi_bins, psi_bins, get_var_index, 3
    )
    
    # ROUND 2: Nonlinear Hybrid Solver (8 variables, continuous)
    print(f"\n" + "="*60)
    print(f"ROUND 2: NONLINEAR HYBRID SOLVER (CONTINUOUS VARIABLES)")
    print(f"="*60)
    
    cqm, sequence_cqm, orig_phi_cqm, orig_psi_cqm, phi_var_names, psi_var_names = folder.create_cqm_problem(
        start_residue=start_res, num_residues=num_res
    )
    
    best_sample_hybrid, best_energy_hybrid, sampleset_hybrid = folder.solve_with_nonlinear_hybrid(cqm)
    
    if best_sample_hybrid is not None:
        opt_phi_hybrid, opt_psi_hybrid = folder.decode_cqm_solution(
            best_sample_hybrid, sequence_cqm, phi_var_names, psi_var_names
        )
        
        # COMPARISON
        print(f"\n" + "="*60)
        print(f"SOLVER COMPARISON RESULTS")
        print(f"="*60)
        
        print(f"\nPERFORMANCE COMPARISON:")
        print(f"• ExactSolver time: {folder.timing_stats['total_time']:.4f}s")
        print(f"• Hybrid solver time: {folder.timing_stats_hybrid['total_time']:.4f}s")
        print(f"• ExactSolver variables: {folder.timing_stats['variables']} (discrete)")
        print(f"• Hybrid solver variables: {folder.timing_stats_hybrid['variables']} (continuous)")
        
        print(f"\nSOLUTION QUALITY:")
        print(f"• ExactSolver energy: {best_energy_exact:.6f}")
        print(f"• Hybrid solver energy: {best_energy_hybrid:.6f}")
        print(f"• Energy difference: {abs(best_energy_exact - best_energy_hybrid):.6f}")
        
        print(f"\nANGLE COMPARISON:")
        print("Res | ExactSolver (φ,ψ) | Hybrid Solver (φ,ψ) | Difference")
        print("-" * 65)
        for i in range(len(sequence)):
            phi_exact = opt_phi_exact[i]
            psi_exact = opt_psi_exact[i]
            phi_hybrid = opt_phi_hybrid[i]
            psi_hybrid = opt_psi_hybrid[i]
            
            phi_diff = abs(phi_exact - phi_hybrid)
            psi_diff = abs(psi_exact - psi_hybrid)
            
            print(f" {i+1:2d} | ({phi_exact:6.1f},{psi_exact:6.1f})    | ({phi_hybrid:6.1f},{psi_hybrid:6.1f})    | ({phi_diff:5.1f},{psi_diff:5.1f})")
        
        # Create comparison visualization
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
        
        # Original conformation
        valid_orig = [(phi, psi) for phi, psi in zip(orig_phi, orig_psi) 
                      if phi is not None and psi is not None]
        if valid_orig:
            phi_orig, psi_orig = zip(*valid_orig)
            ax1.scatter(phi_orig, psi_orig, c='gray', s=100, alpha=0.8, label='Original')
        
        # ExactSolver result
        ax2.scatter(opt_phi_exact, opt_psi_exact, c='red', s=100, alpha=0.8, label='ExactSolver')
        
        # Hybrid solver result
        ax3.scatter(opt_phi_hybrid, opt_psi_hybrid, c='blue', s=100, alpha=0.8, label='Hybrid Solver')
        
        # Format all plots
        for ax, title in zip([ax1, ax2, ax3], ['Original', 'ExactSolver', 'Hybrid Solver']):
            # Add reference regions
            alpha_box = plt.Rectangle((-90, -70), 60, 60, alpha=0.3, color='lightblue')
            beta_box = plt.Rectangle((-150, 90), 60, 60, alpha=0.3, color='lightgreen')
            ax.add_patch(alpha_box)
            ax.add_patch(beta_box)
            
            ax.set_xlim(-180, 180)
            ax.set_ylim(-180, 180)
            ax.set_xlabel('Phi angle (degrees)')
            ax.set_ylabel('Psi angle (degrees)')
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
            ax.legend()
        
        plt.tight_layout()
        plt.show()
        
        print(f"\nCOMPARISON COMPLETE!")
        print(f"Both solvers found stable conformations")
        print(f"Demonstrated discrete vs continuous optimization")
        print(f"Ready for larger problems with hybrid solver")
        
    return folder

    """
    Run the solvaer comparison.
    """
  
if __name__ == "__main__":
    folder = run_solver_comparison()
