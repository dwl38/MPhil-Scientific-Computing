# Written Assignment 1: "A foundation model for atomistic materials chemistry"

The goal of this written assignment was to "explore the basics of machine learning as applied to computational chemistry and materials science, focusing on machine learning interatomic potentials (MLIPs); [...] specifically, to work with a MACE foundation model, which makes detailed simulations accessible without the need for extensive and time-consuming model development. [...] In particular, the student is to apply the model to a system of their choice to understand how it models atomic interactions, and what its strengths and limitations are, [...] to gain practical experience with MLIPs and develop essential skills for applying machine learning in materials science research."

This written assignment was mostly based on the following paper:

- I. Batatia et al. (2023). A foundation model for atomistic materials chemistry. arXiv preprint, 2401.00096. 

with the original intent being to recreate a showcase example from the paper using [MACE-MP-0](https://github.com/ACEsuit/mace-mp/tree/main). However, students were also given the freedom to select systems outside of the examples in the paper, and to use other foundation models if they so desired.


---

In my case, I have chosen to study the thermal decomposition of 1,3,5-triamino-2,4,6-trinitrobenzine (TATB), which is a high explosive with many industrial applications. In particular, when heated under ambient pressure, TATB is known to thermally decompose in a benign manner at roughly ~330&deg;C (i.e. without exploding), hence it is the preferred explosive for applications where safety against unintended detonation is required.

I chose to study this system due to its relatively small molecular size &ndash; thus keeping computational costs low &ndash; and also for its challenging nature against conventional force-field based simulations (e.g. requiring a good model of chemical reactivity amongst organics).

Instead of MACE-MP-0, I chose to use the [MACE-OFF24 'medium' model](https://github.com/ACEsuit/mace-off), which is described in a different paper:

- D.P. Kovács et al. (2023). MACE-OFF23: Transferable Machine Learning Force Fields for Organic Molecules. arXiv preprint, 2312.15211.

due to the better suitability of MACE-OFF in describing organic reactions. The molecular dynamics simulations were performed using [ASE v3.24.0](https://wiki.fysik.dtu.dk/ase/) interfaced with the [MACE calculator](https://github.com/ACEsuit/mace), running either Langevin dynamics (for NVT) or Nosé-Hoover-Parinello-Rahman dynamics (for NpT). 


---

The project structure is as follows:

- `📁 00-common`: Common assets for all following subprojects

    - `📁 geoms`: Starting molecular geometries

    - `📁 mace-models`: Collection of MACE foundation models
    
    - `📁 venv`: Python virtual environment needed to run the following simulations (by executing the respective `run.py` scripts)
 
      | ℹ️ Important |
      |---------------|
      | This repository does not contain the Python environment files themselves, due to large filesizes! Run the script `venv/setup.sh` to download and install the environment, and `source venv/bin/activate` to activate the environment, before running any simulations. |

    - `📜 view.tcl`: Tcl script for VMD to set up view settings (of my personal preference); use as `vmd <file> -e /path/to/view.tcl`

    - `📜 analyze.py`: Python script for analyzing the chemical content of each frame in a trajectory, to tell if and when the TATB molecule(s) have decomposed

- `📁 01-single-molecule`:  A series of NVT simulations of a single TATB molecule at various temperatures

- `📁 02-three-molecules`: A short series of NVT simulations of three TATB molecules in dilute gas phase (abandoned halfway)

- `📁 03-crystal-rtp`: A NpT simulation of a solid TATB crystal at 300K, 1 atm for equilibration
  
- `📁 04-crystal-heat`: A series of NpT simulations of a solid TATB crystal at 1 atm, with increasing temperature
  
- `📁 05-gas-npt`: A series of NpT simulations of what was initially a TATB crystal, but with the wrong boundary conditions, at 1 atm with increasing temperature; intended to represent a "high-density gas"
  
- `📁 06-reaction-neb`: Using NEB to find the energy barrier of the thermal decomposition reaction

> [!NOTE]
> The trajectory files are not in this repository, due to their large sizes; only the simulation logs and processed results are.
