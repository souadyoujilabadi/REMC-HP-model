"""Implementation of the replica exchange Monte Carlo algorithm for protein folding in the HP model.

This program performs a replica exchange Monte Carlo simulation for protein folding in the HP model.
It outputs the lowest current energies and current temperatures of each replica at each REMC step to a file in the results directory, the execution time,
as well as a visualization of the lowest energy conformation found across all replicas and REMC steps.

Usage:
python3 REMC.py -fasta <fasta_file> -numLocalSteps <num_local_steps> -Tmin <min_temp> -Tmax <max_temp> -numReps <num_reps> -maxSteps <max_steps> -is_random <is_random> -seed <seed> -pmWeight <pm_weight>
All the arguments have default values and are optional, except for the fasta file.
"""

# Standard modules
import argparse
import random
import time

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


def main():
    """REMC algorithm for protein folding in the HP model."""
    # Start timer
    start_time = time.time()  

    # Parse command-line arguments
    args = parse_arguments()

    # Seed the random number generator if a seed is provided
    if args.seed is not None:
        random.seed(args.seed)

    # Read the protein sequence from the FASTA file
    protein_sequence = tl.read_fasta(args.fastaFile)

    # Create the initial conformations of the replicas only once before entering the loop (before MC search)
    conformations = tl.initialize_conformations(protein_sequence, args.numReps, args.Tmin, args.Tmax, args.is_random)  # list of Conformation objects

    # Initialize with a high value to capture the lowest
    absolute_lowest_energy = float('inf')  # Initialize to infinity
    best_conformation = None  # Placeholder to store the best conformation

    # Output file to store the lowest energy and temperature of each replica at each REMC step
    with open('results/energies_output.txt', 'w') as outfile:

        flag = 0 # Flag to swap temperatures between replicas

        # REMC loop
        for step in tqdm(range(args.maxSteps), desc="Running REMC"):

            # Perform MC search on existing conformations and swap temperatures between replicas
            for i in range(args.numReps):
                conformations[i] = tl.monte_carlo(args.numLocalSteps, conformations[i], args.pmWeight)  # Update the conformation object in the list of conformations after MC search

            # Write the energy and temperature for each replica (i. e. conformation) to the output file BEFORE swapping temperatures
            for i, conformation in enumerate(conformations):
                outfile.write(f"REMC step {step+1}, Replica {i}: Current energy and temperature: {conformation.energy} {conformation.temperature}\n")

            # Swap temperatures between replicas
            tl.swap_temperatures(conformations, flag)
            # Change the flag value
            flag = 1 if flag == 0 else 0

            # After swapping temperatures between replicas num_replicas times, we run the MC search again on the same replicas (num_replicas conformations)
            # and then swap temperatures again, and so on until the max_steps is reached
            # The loop will stop when the max_steps is reached

            # After swapping temperatures and MC search, identify and print the current lowest energy
            # Update the absolute lowest energy found and store the corresponding conformation
            for conformation in conformations:
                if conformation.energy <= absolute_lowest_energy:
                    absolute_lowest_energy = conformation.energy
                    best_conformation = conformation

        # Compute the execution time
        execution_time = time.time() - start_time
        print(f"Execution time: {execution_time:.2f} seconds")
        # Print the absolute lowest energy
        print(f"Lowest energy across all replicas and REMC steps: {absolute_lowest_energy}")  
        outfile.write(f"\nExecution time: {execution_time:.2f} seconds")

    # Visualize the best conformation across all replicas and REMC steps
    if best_conformation:
        tl.visualize_2d_conformation(best_conformation)

if __name__ == "__main__":  # If the program REMC.py is executed as a script in a shell, the result of the if test will be True and the corresponding instruction block will be executed.
    main()                  # Otherwise, if the program REMC.py is imported as a module, then the result of the test if will be False (and the corresponding instruction block will not be executed).
