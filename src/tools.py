"""This module contains functions for performing replica exchange Monte Carlo
simulations in the HP model, as well as visualization of the simulations.

Classes:
    Conformation
"""

# Standard modules
import copy
import random
import math
import matplotlib.pyplot as plt
import numpy as np

# Local modules
from conformation import Conformation


def read_fasta(filename):
    """Read a protein sequence from a FASTA file.

    Args:
        filename (str): The path to the FASTA file.

    Returns:
        str: The protein sequence.
    """
    protein_sequence = ""
    with open(filename, "r") as file:
        for line in file:
            if line.startswith(">"):
                continue  # Skip the header line
            protein_sequence += line.strip()
    return protein_sequence


def monte_carlo(numLocalSteps, conformation, pm_weight):
    """Perform a Monte Carlo simulation for conformational sampling.

    This function implements a Monte Carlo simulation to explore different conformations.
    It starts from an initial conformation and performs 'numLocalSteps' steps,
    selecting movements based on 'pm_weight'.

    Parameters:
        - numLocalSteps (int): The number of local steps to perform in the Monte Carlo simulation.
        - conformation (Conformation): The initial conformation to start the simulation from.
        - pm_weight (float): Control parameter for movement selection.

    Returns:
    - Conformation: The final conformation after the Monte Carlo simulation.
"""

    K_b = 0.0019872  # Boltzmann constant in kcal/mol*K
    sequence_length = len(conformation.sequence)

    for _ in range(numLocalSteps):
        # Create a copy of the current conformation 
        new_conformation = copy.deepcopy(conformation)  # Deep copy, so the new conformation is a different object from the old one

        # Choose a random residue to perform a movement on
        residue_index = random.randint(0, sequence_length - 1)

        # Choose a movement for the selected residue based on the pm_weight
        new_conformation.choose_movement_randomly(residue_index, pm_weight)

        # Compute the energy of the new and current conformations
        new_energy = new_conformation.compute_energy()
        energy = conformation.compute_energy()

        # Metropolis criterion to decide whether to accept the new conformation
        proba = math.exp((-(new_energy - energy))/(conformation.temperature * K_b))

        if new_energy < energy or random.uniform(0,1) <= proba:
            conformation = new_conformation

    return conformation


def initialize_conformations(protein_sequence, num_replicas, min_temperature, max_temperature, is_random):
    """Initialize conformations for replica exchange simulations.

    This function initializes a list of conformations, in which each conformation belongs to a replica (i. e. represents a replica).
    Each replica corresponds to a temperature between 'min_temperature' and 'max_temperature'.
    Initial conformations can be random if 'is_random' is True, otherwise, it's linear.

    Parameters:
        protein_sequence (str): Protein sequence.
        num_replicas (int): Number of replicas.
        min_temperature (float): Minimum temperature.
        max_temperature (float): Maximum temperature.
        is_random (bool): Assign random initial conformations if True.

    Returns:
        conformations (list): List of Conformation (initialized conformations).
    """
    conformations = []  # List to store initialized conformations

    # Calculate the temperature step size for evenly distributed temperatures
    temperature_step = (max_temperature - min_temperature) / (num_replicas - 1) 

    for i in range(num_replicas):

        # Calculate the temperature for the current replica
        temperature = min_temperature + i * temperature_step

        # Create a new conformation instance for the current temperature
        conformation = Conformation(protein_sequence, temperature)

        # Assign an initial conformation, either randomly or linearly
        conformation.assign_initial_conformation(is_random)

        # Add the initialized conformation to the list
        conformations.append(conformation)

    return conformations  # list of Conformation objects


def swap_temperatures(conformations, flag):
    """Attempt temperature swaps between adjacent replicas in a replica exchange simulation.

    This function attempts to swap temperatures between adjacent replicas in a replica exchange simulation.
    The probability of swapping is determined by the Metropolis criterion based on the temperature
    difference and energy difference between replicas.

    Parameters:
        conformations (list of Conformation): List of conformations representing different replicas.
        flag (int): Flag indicating the starting index for swapping attempts.
    """

    # Get the total number of replicas
    n_replica = len(conformations)

    # Define a scaling factor for temperature calculations
    K_b2 = 0.0001679010

    # Initialize the index for swapping attempts
    i = flag

    # Loop through the replicas for swapping attempts
    while i < (n_replica - 1):
        j = i + 1

        # Calculate inverse temperatures and energies for the two replicas
        betai = 1 / (K_b2 * conformations[i].temperature)
        betaj = 1 / (K_b2 * conformations[j].temperature)
        energyi = conformations[i].compute_energy()
        energyj = conformations[j].compute_energy()

        # Calculate the energy and temperature difference for the Metropolis criterion
        delta = (betaj - betai) * (energyi - energyj)

        # Calculate the probability of swapping temperatures
        proba = math.exp(-delta)  # Metropolis criterion to decide if the temperatures should be swapped

        if delta > 0 or random.uniform(0, 1) < proba:
            # Swap the temperatures of the adjacent replicas
            conformations[i].temperature, conformations[j].temperature = conformations[j].temperature, conformations[i].temperature

        # Increment the index to consider the next pair of adjacent replicas
        i += 2


def visualize_2d_conformation(conformation):
    """Generate a 2D representation of a conformation.

    Parameters:
        conformation (Conformation): The conformation to visualize.
    """

    # Get the positions of the residues in the conformation
    positions = conformation.positions

    # Get the x and y coordinates of the positions
    x_coords = [pos[0] for pos in positions]
    y_coords = [pos[1] for pos in positions]

    # Split the positions based on the HP model
    hydrophobe_positions = [pos for i, pos in enumerate(positions) if conformation.residues[i].polarity_classification == 'H']
    polar_positions = [pos for i, pos in enumerate(positions) if conformation.residues[i].polarity_classification != 'H']

    # Plot the hydrophobe and polar points with labels for legend
    plt.scatter([pos[0] for pos in hydrophobe_positions], [pos[1] for pos in hydrophobe_positions], color='black', s=100, label='Hydrophobe')
    plt.scatter([pos[0] for pos in polar_positions], [pos[1] for pos in polar_positions], color='red', s=100, label='Polar')

    # Plot the lines to connect the points
    plt.plot(x_coords, y_coords, color='blue')

    # Label the points based on the HP model
    for i, (x, y) in enumerate(positions):
        label = f"{i}"
        plt.annotate(label, (x, y), color='white', ha='center', va='center', fontsize=8)

    # Set x and y ticks based on the min and max values of the coordinates to align the grid
    plt.xticks(np.arange(min(x_coords)-1, max(x_coords)+2, 1.0))  # start, stop, step
    plt.yticks(np.arange(min(y_coords)-1, max(y_coords)+2, 1.0))

    plt.title("2D Visualization of the Lowest Energy Conformation (HP model)")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.savefig('results/2d_conformation.png')
