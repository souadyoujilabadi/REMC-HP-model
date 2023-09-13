"""Implementation of the replica exchange Monte Carlo algorithm for protein folding in the HP model.

This program performs a replica exchange Monte Carlo simulation for protein folding in the HP model.
The simulation is parallelized using the multiprocessing module to speed up the computation.
The parallelization is done at the level of the Monte Carlo search, i.e. each replica is simulated in parallel.
It outputs the lowest current energies and current temperatures of each replica at each REMC step to a file in the results directory, the execution time,
as well as a visualization of the lowest energy conformation found across all replicas and REMC steps.

Usage:
python3 REMC_parallel.py -fasta <fasta_file> -numLocalSteps <num_local_steps> -Tmin <min_temp> -Tmax <max_temp> -numReps <num_reps> -maxSteps <max_steps> -is_random <is_random> -seed <seed> -pmWeight <pm_weight>
All the arguments have default values and are optional, except for the fasta file.
"""

# Standard modules
import argparse
import random
import time
from multiprocessing import Pool

# Local modules
import tools as tl
from tqdm import tqdm


def parse_arguments():
    """Parse command-line arguments.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments and their values.
    """
    parser = argparse.ArgumentParser(description="2D protein folding simulation using REMC method with Dill's HP model")
    parser.add_argument("-fasta", type=str, dest="fastaFile", help="Protein fasta file", required=True)
    parser.add_argument("-numLocalSteps", type=int, dest="numLocalSteps", default=500, help="Number of local steps in a Monte Carlo search")
    parser.add_argument("-Tmin", type=float, dest="Tmin", default=160.0, help="Minimum temperature value")
    parser.add_argument("-Tmax", type=float, dest="Tmax", default=220.0, help="Maximum temperature value")
    parser.add_argument("-numReps", type=int, dest="numReps", default=5, help="Number of replicas to simulate")
    parser.add_argument("-maxSteps", type=int, dest="maxSteps", default=10, help="Maximum number of steps")
    parser.add_argument("-is_random", type=str_to_bool, dest="is_random", default=True, help="Whether to assign a random initial conformation (True) or linear (False)")
    parser.add_argument("-seed", type=int, dest="seed", default=None, help="Seed for random number generation")
    parser.add_argument("-pmWeight", type=float, dest="pmWeight", default=0.4, help="Weight to give pull moves vs. VSHD moves")
    return parser.parse_args()


def str_to_bool(value):
    """Convert a string to a boolean."""
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# This is the wrapper function
# It is needed because the multiprocessing module requires a top-level function
def parallel_monte_carlo(args):
    """Wrapper function for the monte_carlo function."""
    # Unpack the arguments and call the monte_carlo function
    return tl.monte_carlo(*args)  # *args is used to unpack the arguments


def main():
    """REMC algorithm for protein folding in the HP model using multiprocessing."""
    # Start timer
    start_time = time.time()  

    # Parse command-line arguments
    args = parse_arguments()

    # Seed the random number generator if a seed is provided
    if args.seed is not None:
        random.seed(args.seed)

    # Read the protein sequence from the FASTA file
    protein_sequence = tl.read_fasta(args.fastaFile)

    # Create the initial conformations
    conformations = tl.initialize_conformations(protein_sequence, args.numReps, args.Tmin, args.Tmax, args.is_random)

    # Placeholders for the lowest energy and best conformation found
    absolute_lowest_energy = float('inf')  # Initialize to infinity
    best_conformation = None

    # Output file
    with open('results/energies_output.txt', 'w') as outfile:

        # Use multiprocessing for parallel execution on the replicas
        with Pool(processes=args.numReps) as pool:
            flag = 0  # Flag to swap temperatures between replicas

            # Run the REMC algorithm for the specified number of steps
            for step in tqdm(range(args.maxSteps), desc="Running REMC"):
                # Parallelize the MC searches across all replicas
                conformations = pool.map(parallel_monte_carlo, [(args.numLocalSteps, conformation, args.pmWeight) for conformation in conformations])

                # Write the current energies and temperatures of all replicas to the output file
                for i, conformation in enumerate(conformations):
                    outfile.write(f"REMC step {step+1}, Replica {i}: Current energy and temperature: {conformation.energy} {conformation.temperature}\n")

                # Swap temperatures between replicas
                tl.swap_temperatures(conformations, flag)
                # Change the flag value
                flag = 1 if flag == 0 else 0

                # Update the lowest energy and best conformation found
                for conformation in conformations:
                    # Check if the current conformation has the lowest energy
                    if conformation.energy <= absolute_lowest_energy:
                        # Update the lowest energy and best conformation found
                        absolute_lowest_energy = conformation.energy
                        best_conformation = conformation

        # Stop timer
        execution_time = time.time() - start_time
        print(f"Execution time: {execution_time:.2f} seconds")
        print(f"Lowest energy across all replicas and REMC steps: {absolute_lowest_energy}")

        # Write the execution time to the output file
        outfile.write(f"\nExecution time: {execution_time:.2f} seconds")


    # Visualize the best conformation found
    if best_conformation:
        tl.visualize_2d_conformation(best_conformation)

if __name__ == '__main__':
    main()
