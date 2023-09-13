# REMC-HP-model

This is a Python implementation of the Replica Exchange Monte Carlo (REMC) algorithm for protein folding in the Hydrophobic-Polar (HP) model. The program uses a parallelized version of the REMC algorithm to speed up the computation.

## Requirements

To run the program, you need to setup your environment:

- Clone the repository:
  
```bash
git clone https://github.com/souadyoujilabadi/REMC-HP-model.git
```

- Move to the new directory:

```bash
cd REMC-HP-model
```

- Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)

- Create the `env_REMC` conda environment:
  
```bash
conda env create -f env_REMC.yml
```

- Load the `env_REMC` conda environment:
  
```bash
conda activate env_REMC
```

## Usage

To run the program, you need to specify the input parameters using command-line arguments. Here's an example usage:
```bash
python3 REMC.py -fasta input.fasta 
```

```bash
python3 REMC_parallel.py -fasta input.fasta
```

A sequence was generated to be equivalent to the sequence tested by Chris Thachuk et al. in the provided article (see README in the article's provided code), and used in order to compare the results of the REMC algorithm implemented in this project with the results of the REMC algorithm implemented by Chris Thachuk et al:

```bash
python3 REMC.py -fasta test.fasta.txt
```

```bash
python3 REMC_parallel.py -fasta test.fasta.txt
```

The available command-line arguments are:

|Parameter      | Description                                                              | Default value |
|_______________|__________________________________________________________________________|_______________|
|-fasta         | The path to the protein sequence FASTA file.                             |               |
|-numLocalSteps | The number of local steps in a Monte Carlo search.                       | 500           |
|-Tmin          | The minimum temperature value.                                           | 160.0         |
|-Tmax          | The maximum temperature value.                                           | 220.0         |
|-numReps       | The number of replicas to simulate.                                      | 5             |
|-maxSteps      | The maximum number of steps.                                             | 10            |
|-is_random     | Whether to assign a random initial conformation (True) or linear (False).| True          |
|-seed          | Seed for random number generation (results reproducibility).             | None          |
|-pmWeight      | Weight to give pull moves vs. VSHD moves.                                | 0.4           |

All the arguments have default values and are optional, except for the fasta file.

## Output

The program outputs a file `energies_output.txt` that contains the lowest current energy for each replica at the current temperature.
The program also generates a plot of the conformation with the lowest energy across all replicas and REMC steps.
