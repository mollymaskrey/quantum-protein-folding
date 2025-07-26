
# Quantum vs Classical Protein Folding Optimization

A comparative study demonstrating where quantum computers become essential for protein folding optimization problems that are intractable for classical hardware.

## 🧬 Overview

This repository contains the complete experimental pipeline comparing classical exhaustive search (ExactSolver) with quantum annealing (D-Wave Hybrid) for protein dihedral angle optimization. Our study identifies the computational transition point where quantum advantage becomes necessary.

### Key Findings
- **24 variables**: Classical solver completes in ~9 seconds, quantum in ~1.8 seconds
- **28 variables**: Classical solver fails after 15+ minutes, quantum completes in ~1.4 seconds
- **Transition point**: Quantum computers become essential between 24-28 variables
- **Hardware tested**: MacBook Pro M4 Max (128GB RAM) vs D-Wave quantum cloud

## 🚀 Quick Start

### Prerequisites
- Python 3.9 or higher
- D-Wave account and API token (for quantum experiments)
- 8GB+ RAM recommended for classical experiments

### Installation

1. Clone this repository:
```bash
git clone https://github.com/[your-username]/quantum-protein-folding.git
cd quantum-protein-folding
```

2. Create and activate a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

4. Set up D-Wave credentials (for quantum experiments):
```bash
export DWAVE_SAMPLER_TOKEN="your-dwave-token-here"
```

### Running the Experiments

#### 1. Main Comparison Study (24 variables)
Compare classical and quantum approaches on the same 4-residue protein segment:

```bash
python src/quantum_protein_comparison.py
```

**Expected output:**
- Downloads human insulin structure from AlphaFold
- Extracts dihedral angles and protein sequence
- Runs classical ExactSolver (~9 seconds)
- Runs quantum hybrid solver (~1.8 seconds)
- Generates comparative analysis and visualizations

#### 2. Quantum Advantage Test (28 variables)
Demonstrate quantum solving of classically intractable problems:

```bash
python src/quantum_28_variable_test.py
```

**Expected output:**
- Tests 7-residue segment with backbone and side chain angles
- Shows quantum solution where classical methods fail
- Demonstrates consistent quantum performance scaling

## 📊 Results

### Performance Comparison

| Problem Size | Classical Time | Quantum Time | Speedup |
|--------------|----------------|--------------|---------|
| 24 variables | 8.99 seconds   | 1.77 seconds | 5.1x    |
| 28 variables | Failed (>15min)| 1.40 seconds | ∞       |

### Hardware Specifications

**Classical Computing:**
- MacBook Pro with M4 Max chip
- 128 GB RAM
- Local computation

**Quantum Computing:**
- D-Wave quantum annealing system
- Cloud-based hybrid solver
- Continuous variable optimization

## 🧪 Scientific Impact

This study demonstrates:

1. **Classical computational limits**: High-end hardware reaches practical limits at 28 variables
2. **Quantum consistency**: Quantum performance remains stable across problem sizes
3. **Transition point identification**: Clear boundary where quantum becomes necessary
4. **Real protein applications**: Using actual biological data (human insulin)

## 📁 Repository Structure

```
quantum-protein-folding/
├── README.md                           # This file
├── requirements.txt                    # Python dependencies
├── solver_comparison.py
├── run_28_variable_quantum_test.py
├── docs/
│    └── paper.pdf                       # Full research paper
└── results/
    ├── solver_comparison               # Results for up to 24 Variables
    └── run_28_variable_quantum_test    # 28 Variable Results using nonlinear hybrid solver
```

## 🔬 Methodology

### Classical Approach (ExactSolver)
- **Method**: Exhaustive search with discrete variables
- **Variables**: Binary representation of dihedral angle bins
- **Scaling**: Exponential (2^n solutions)
- **Advantage**: Finds guaranteed global optimum within discrete space

### Quantum Approach (D-Wave Hybrid)
- **Method**: Quantum annealing with continuous variables
- **Variables**: Real-valued dihedral angles (-180° to +180°)
- **Scaling**: Consistent performance across problem sizes
- **Advantage**: Handles problems impossible for classical computers

### Energy Function
Both methods optimize protein conformational energy based on:
- Ramachandran plot preferences (alpha helix, beta sheet regions)
- Steric clash avoidance
- Amino acid-specific constraints (e.g., proline restrictions)

## 📈 Usage Examples

### Basic Protein Optimization
```python
from src.quantum_protein_comparison import QuantumProteinFolder

# Initialize the solver
folder = QuantumProteinFolder()

# Download and analyze a protein
pdb_file = folder.download_alphafold_protein('P01308')  # Human insulin
folder.analyze_protein_structure(pdb_file)
folder.extract_dihedral_angles()

# Run quantum optimization
cqm, sequence, phi_vars, psi_vars = folder.create_cqm_problem(start_residue=10, num_residues=4)
sample, energy, sampleset = folder.solve_with_nonlinear_hybrid(cqm)
```

### Custom Protein Analysis
```python
# Use your own protein structure
folder = QuantumProteinFolder()
folder.analyze_protein_structure('my_protein.pdb')

# Optimize specific regions
results = folder.optimize_segment(start=20, length=6)
```

## Configuration

### D-Wave Setup
1. Create account at [D-Wave Leap](https://cloud.dwavesys.com/leap/)
2. Get your API token from the dashboard
3. Set environment variable or configure in code:
```python
import os
os.environ['DWAVE_SAMPLER_TOKEN'] = 'your-token-here'
```

### Computational Requirements
- **Minimum**: 128GB RAM, Python 3.8+
- **Quantum access**: D-Wave Leap account 

## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@article{yourname2024quantum,
  title={Quantum Computing Meets Protein Folding: A Side-by-Side Comparison of Classical and Quantum Optimization Methods},
  author={Your Name},
  journal={Journal Name},
  year={2024},
  note={Code available at: https://github.com/[your-username]/quantum-protein-folding}
}
```

## 🤝 Contributing

We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) for details.

### Areas for Contribution
- Additional protein test cases
- Alternative quantum solvers (IBM Qiskit, etc.)
- Enhanced energy functions
- Performance optimizations
- Documentation improvements

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🔗 Related Work

- [D-Wave Ocean SDK](https://github.com/dwavesystems/dwave-ocean-sdk)
- [AlphaFold Protein Database](https://alphafold.ebi.ac.uk/)
- [BioPython](https://biopython.org/)
- [Quantum Optimization Research](https://docs.dwavesys.com/docs/latest/handbook_applications.html)

## 📞 Contact

- **Author**: Your Name
- **Email**: your.email@domain.com
- **Institution**: Your University/Organization
- **ORCID**: [0000-0000-0000-0000](https://orcid.org/0000-0000-0000-0000)

## 🙏 Acknowledgments

- D-Wave Systems for quantum computing access
- AlphaFold team for protein structure data
- BioPython developers for structural analysis tools
- Quantum computing research community

---

**Keywords**: quantum computing, protein folding, optimization, D-Wave, quantum annealing, computational chemistry, drug discovery
