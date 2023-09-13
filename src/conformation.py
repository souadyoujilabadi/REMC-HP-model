"""A module that represents protein conformations based on a sequence of residues.

This module provides the `Conformation` class which models a protein conformation
in a 2D spatial arrangement.

Classes:
    Conformation: Represents a protein conformation and provides methods to manage its spatial arrangement and to compute its energy.
"""

# Standard modules
import random
import math

# Local modules
from residue import Residue


class Conformation:
    """Represents a protein conformation.

    Attributes:
        sequence (str): The protein's sequence as a string of residue characters.
        residues (list): List of Residue objects representing the residues in the sequence.
        positions (list of lists): 2D coordinates representing each residue's spatial location.
        temperature (float): The temperature of the protein conformation.
        energy (float): The energy of the protein conformation.
    """

    def __init__(self, sequence, temperature):
        """Initialize a Conformation object with a given sequence and temperature.

        Parameters:
            sequence (str): Sequence of residues.
            temperature (float): Operational temperature for the protein conformation.
        """
        self.sequence = sequence
        self.residues = self.sequence_to_residue_objects()  # list of Residue objects
        self.positions = self.initialize_positions()  # list of 2D coordinates for each residue in the sequence
        self.temperature = temperature
        self.energy = 0  # default energy value to be updated after energy computation

    def sequence_to_residue_objects(self):
        """Convert the protein sequence into a list of Residue objects.
        The method iterates over the sequence and creates a Residue object for each residue.

        Returns:
            list: List of Residue objects that maintain the same order as the residues in the sequence.
        """
        residues = []
        for residue in self.sequence:
            residue_obj = Residue(residue)  # Instantiate class Residue = Create Residue object
            residues.append(residue_obj)  # Append the Residue object to the residues list
        return residues

    def initialize_positions(self):
        """Initialize the positions list with default 2D coordinates.
        This method is primarily used to reset the positions of residues during random placement when no adjacent position is found.

        Returns:
            list of lists: List of 2D coordinates [[x, y], [x, y], ...] for each residue in the sequence.
        """
        positions = [[0, 0] for position in self.sequence]
        return positions

    def find_adjacent_empty_position(self, residue_index):
        """Find an unoccupied position adjacent to the specified residue.

        Parameters:
            residue_index (int): Index of the given residue.

        Returns:
            list: List containing the available coordinates [x, y], if a free position is found. Otherwise, it returns nothing.
        """
        # Define possible adjacent moves relative to a given position.
        possible_moves = [[-1, 0], [1, 0], [0, -1], [0, 1]]
        random.shuffle(possible_moves)  # Randomize the order in which we check the moves

        # Get the current position of the residue at the given index
        current_position = self.positions[residue_index]

        for move in possible_moves:
            # Get the adjacent position based on the current position and the move
            adj_position = [current_position[0] + move[0], current_position[1] + move[1]]  
            if adj_position not in self.positions:
                return adj_position
        return None

    def find_diagonal_empty_position(self, residue_index1, residue_index2):
        """Find an unoccupied diagonal positions relative to the specified residues.

        Parameters:
            residue_index1 (int): Index of the first residue in the sequence.
            residue_index2 (int): Index of the second residue in the sequence.

        Returns:
            tuple: A tuple containing two lists, each representing the unoccupied diagonal
            positions relative to residue_index1 and residue_index2, respectively.
            If no such positions are found, it return nothing.
        """
        # Define possible diagonal moves relative to a given position.
        possible_pull_moves = [[-1, -1], [-1, 1], [1, -1], [1, 1]]
        random.shuffle(possible_pull_moves)  # Shuffle to ensure random selection of moves

        # Get the current positions of the specified residues
        current_position1 = self.positions[residue_index1]
        current_position2 = self.positions[residue_index2]

        # Iterate through possible diagonal moves to find unoccupied positions
        for move in possible_pull_moves:
            # Calculate the new diagonal positions relative to the current positions of the residues
            adj_position1 = [current_position1[0] + move[0], current_position1[1] + move[1]]
            adj_position2 = [current_position2[0] + move[0], current_position2[1] + move[1]]

            # Check if the calculated positions are unoccupied
            if adj_position1 not in self.positions and adj_position2 not in self.positions:
                return adj_position1, adj_position2
        return None, None

    def move(self, adj_position, residue_index):
        """Update the position of the specified residue.

        Parameters:
            adj_position (list): List containing unoccupied coordinates [x, y].
            residue_index (int): Index of the residue in the sequence.

        Returns:
            None
        """
        self.positions[residue_index] = adj_position

    def assign_initial_conformation(self, is_random=True):
        """Assign initial 2D coordinates to residues based on a desired arrangement.
        This method provides two possible arrangements for setting the spatial coordinates of protein residues:

        1. Linear Arrangement: Residues are aligned in a straight horizontal line, without any twists or turns.
        2. Random Arrangement: Each residue is placed in an available position adjacent to the previous one.
        If no available spot is found, the method backtracks to the previous residue to find an alternative position. 
        If it reaches the first residue (Nt) and no available position is found, it resets all positions and retries.

        Parameters:
            is_random (bool, default=True): Determines the arrangement type:
                                            If True, residues are placed randomly. Otherwise, they're placed linearly.
                                            
        Returns:
            None: The coordinates of residues are updated and stored in the 'positions' attribute of the class.

        NB:
            The infinite 2D grid assumes no explicit boundaries, and the starting point is considered [0, 0].
        """
        sequence_length = len(self.sequence)
        # Define the initial placement position for the first residue
        placement_origin = [0, 0]  # In a 2D grid without explicit boundaries, the grid can be seen as an infinite 2D plane)
        
        if is_random:
            # Random arrangement

            # Position the first residue at the origin
            self.positions[0] = placement_origin 

            # Handle the positioning of the remaining residues
            idx = 1  # Start placing from the second residue 
            while idx < sequence_length:
                # Find an empty position adjacent to the last placed residue
                adj_position = self.find_adjacent_empty_position(idx - 1)

                if adj_position:  # If an adjacent position is found
                    # If a position is found, assign it to the current residue
                    self.positions[idx] = adj_position
                    idx += 1 # Move on to the next residue
                else:
                    # Backtrack to the previous residue and try again
                    if idx > 1:  # Ensure we're not already at the start
                        idx -= 1  # Go back to the previous residue
                    else:
                        # If we're back at the beginning and can't find a spot, start over (=reset all positions)
                        self.positions = self.initialize_positions()
                        idx = 1  # Restart from the second residue

        else:
            # Linear arrangement
            for idx in range(sequence_length):
                # Set the residue's position linearly along the x-axis, updating the list of positions
                self.positions[idx] = [placement_origin[0] + idx, placement_origin[1]]

    def end_move(self, residue_index):
        """Move a residue located at the beginning or the end of the sequence
        to an adjacent available position if one exists.

        This method checks if the residue at the given index is at the start
        or the end of the sequence. If it is, and there's an available position
        next to it, the residue is moved to that position. The success or failure
        of the movement is indicated by the return value.

        Parameters:
            residue_index (int): The index of the residue to be moved.

        Returns:
            bool: True if the movement was successful, False otherwise.
        """
        # Placeholder for potential adjacent position
        adj_position = None

        # If the residue is at the beginning of the sequence
        if residue_index == 0:
            # Attempt to find an empty position next to the second residue
            adj_position = self.find_adjacent_empty_position(residue_index + 1)

        # If the residue is at the end of the sequence
        elif residue_index == len(self.sequence) - 1:
            # Attempt to find an empty position next to the penultimate residue
            adj_position = self.find_adjacent_empty_position(residue_index - 1)

        # If an adjacent position was found
        if adj_position is not None:
            # Relocate the residue to the discovered adjacent position
            self.move(adj_position, residue_index)
            return True  # Movement was successful

        # Return False if no valid movement was made
        return False

    def corner_move(self, residue_index):
        """Attempt a corner move on the specified residue.

        This move shifts a residue to an adjacent corner. The move's feasibility depends on
        the positions of the neighboring residues and whether the destination corner is already occupied.

        Parameters:
            residue_index (int): The target residue's index.

        Returns:
            bool: True if the movement was successful, False otherwise. Failures occur if:
                  1. The residue is not at either end of the sequence.
                  2. The neighboring residues aren't in a configuration allowing a diagonal move.
                  3. The desired corner position is already occupied.
        """

        # Avoid moves for first or last residues
        if residue_index == 0  or residue_index == (len(self.sequence) - 1):
            return False

        prev_residue_coord = self.positions[residue_index - 1]
        current_residue_coord = self.positions[residue_index]
        next_residue_coord = self.positions[residue_index + 1]

        # Ensure the neighboring residues form a configuration that allows a diagonal move
        if math.dist(prev_residue_coord, next_residue_coord) != math.sqrt(2):  # diagonal of a unit square = sqrt(2) (pythagorean theorem)
            return False

        # Calculate the new potential position for the residue
        # The logic is based on the alignment of the current residue with its neighbors.
        # Specifically, if the current and next residues are vertically aligned (same x-coordinate), the residue should move horizontally. 
        # Conversely,if they are horizontally aligned, the residue should move vertically.
        new_position = [
            prev_residue_coord[0] if next_residue_coord[0] == current_residue_coord[0] else next_residue_coord[0],
            prev_residue_coord[1] if prev_residue_coord[1] == current_residue_coord[1] else next_residue_coord[1]
            ]

        # Execute the move only if the target corner is not already occupied
        if new_position not in self.positions:
            self.move(new_position, residue_index)
            return True  # Movement was successful

        return False  # Movement failed

    def crankshaft_move(self, residue_index):
        """Attempts a crankshaft move on the residues at the given index and its subsequent index.

        The crankshaft move involves rotating two residues around their respective adjacent residues:
        - Residue at index `residue_index` pivots around the residue at index `residue_index-1`.
        - Residue at index `residue_index+1` pivots around the residue at index `residue_index+2`.

        This rotation only occurs if:
        1. The chosen residue isn't the first, the second-to-last, or the last residue in the sequence.
        2. The distance between the pivoting residues (i. e. the residues around which the move takes place: `residue_index-1` and `residue_index+2`) is 1.
        3. The new positions resulting from the move are unoccupied.

        Parameters:
            residue_index (int): The index of the residue where the crankshaft move begins.

        Returns:
            bool: True if the movement was successful, False otherwise.
        """

        # Ensure residue_index is not the first, the second-to-last, or the last residue in the sequence
        if residue_index == 0 or residue_index >= (len(self.sequence) - 2):
            return False

        # Determine the pivot points around which the crankshaft rotation will occur
        pivot_start = self.positions[residue_index - 1]
        pivot_end = self.positions[residue_index + 2]

        # Ensure the pivoting points are connected, required for a crankshaft move
        if math.dist(pivot_start, pivot_end) != 1:
            return False

        # Calculate the new coordinates after the crankshaft rotation
        rotated_pos1 = [
            pivot_start[0] + (pivot_start[0] - self.positions[residue_index][0]), # calculates the difference in the x-coordinates between pivot_start and the chosen residue, then adds this difference to the x-coordinate of pivot_start to get the new x-coordinate for the residue at residue_index = move in the opposite direction of the chosen residue
            pivot_start[1] + (pivot_start[1] - self.positions[residue_index][1]) # same for the y-coordinate
        ]
        rotated_pos2 = [
            pivot_end[0] + (pivot_end[0] - self.positions[residue_index + 1][0]), # same considering the subsequent residue (residue_index + 1)
            pivot_end[1] + (pivot_end[1] - self.positions[residue_index + 1][1])
        ]

        # Ensure the calculated positions aren't occupied by other residues
        if rotated_pos1 not in self.positions and rotated_pos2 not in self.positions:
            # Execute the move if both positions are free
            self.move(rotated_pos1, residue_index)
            self.move(rotated_pos2, residue_index + 1)
            return True  # Movement was successful

        return False  # Rotation failed due to occupied positions

    def is_connected(self, residue_index1, residue_index2):
        """Checks if two residues are connected.

        Parameters:
            residue_index1 (int): Index of the first residue.
            residue_index2 (int): Index of the second residue.

        Returns:
            bool: True if the residues are connected; False otherwise."""

        # Relative positions for connected neighbors
        connected_neighbors = [[-1, 0], [1, 0], [0, -1], [0, 1]]

        for neighbor in connected_neighbors:
            # Calculate potential position of a neighboring residue
            potential_position = [self.positions[residue_index1][0] + neighbor[0],
                                  self.positions[residue_index1][1] + neighbor[1]]
            # Return True if positions match
            if potential_position == self.positions[residue_index2]:
                return True

        return False

    def pull_move(self, residue_idx):
        """Attemps a pull move.
        This method adjusts the position of a given residue in a sequence based on a randomly chosen direction.
        It checks if the given residue and its neighboring residue can be pulled in a diagonal direction
        while maintaining the sequence connectivity.

        The residue's position will be updated if it satisfies the following conditions:
        1. The residue isn't located at the ends of the sequence.
        2. The new positions for the current residue and its neighboring residue are not the same.
        3. The distance between the new positions for the residues is 1.
        4. The diagonal distance between the current and the new position for each residue is sqrt(2).

        Parameters:
            residue_idx (int): Index of the residue starting the pull move.

        Returns:
            bool: True if the movement was successful, False otherwise.

        Note: The function updates the positions of residues in place and doesn't return the updated sequence."""

        # Define adjustments for each direction (movement patterns based on the direction)
        pull_directions = {
            'towards_end': {
                'next_residue_index': residue_idx + 1,
                'previous_residue_index': residue_idx - 1,
                'forward_adjustment_step': 2,
                'backward_adjustment_step': 1
            },
            'towards_start': {
                'next_residue_index': residue_idx - 1,
                'previous_residue_index': residue_idx + 1,
                'forward_adjustment_step': -2,
                'backward_adjustment_step': -1
            }
        }

        # Randomly select a direction either 'towards_end' or 'towards_start' for pulling the residue.
        direction = random.choice(list(pull_directions.keys()))  # Either 'towards_end' or 'towards_start'
        adjustments = pull_directions[direction]  # Get the adjustments (movement patterns) for the chosen direction

        # Return False if residue is at either end of the sequence
        if residue_idx == 0 or residue_idx == (len(self.sequence) - 1):
            return False

        # Identify new diagonal positions for the current and next residue based on the chosen direction
        new_pos_for_current, new_pos_for_next = self.find_diagonal_empty_position(residue_idx, adjustments['next_residue_index'])

        # Extract the current positions of the residues in focus
        pos_for_current = self.positions[residue_idx]
        pos_for_next = self.positions[adjustments['next_residue_index']]

        # Check if the new positions satisfy the constraints for pulling in the diagonal direction
        if (new_pos_for_current and new_pos_for_next and
            new_pos_for_current != new_pos_for_next and
            math.dist(new_pos_for_current, new_pos_for_next) == 1 and
            math.dist(new_pos_for_current, pos_for_next) == math.sqrt(2) and
            math.dist(new_pos_for_next, pos_for_current) == math.sqrt(2)):

            # Update positions of the residues based on the calculated new positions
            self.move(new_pos_for_next, adjustments['next_residue_index'])
            self.move(new_pos_for_current, residue_idx)

            # Loop through subsequent residues to ensure the sequence remains connected
            i = adjustments['forward_adjustment_step']
            while 0 <= residue_idx - i < len(self.sequence) and not self.is_connected(residue_idx, residue_idx - i):
                temp_pos = self.positions[residue_idx - i]
                self.move(pos_for_next, residue_idx - i)
                pos_for_next = temp_pos
                i += adjustments['backward_adjustment_step']
                return True

        return False

    def weighted_random_choice(self, choices, weights):
        """Select a random choice based on the provided weights.

        This method returns a randomly-selected element from the 'choices' list,
        where the probability of each element being chosen is determined by the
        corresponding weights in the 'weights' list.

        Parameters:
            choices (list): A list of items to choose from.
            weights (list): A list of weights corresponding to each item in 'choices'. The weights sum up to 1.0.

        Returns:
            str: A randomly-selected item from 'choices'.

        Example:
            >>> weighted_random_choice(['corner', 'crankshaft', 'pull'], [0.3, 0.3, 0.4])
            'pull'  # this is just one possible output given the weights
        """
        return random.choices(choices, weights=weights, k=1)[0]  # k=1 : return a single element from the list of choices. [0] : extract this item from the list.

    def choose_movement_randomly(self, residue_index, pm_weight):
        """Modifies the conformation of the protein based on the provided residue index.

        The type of movement depends on the position of the residue:
        - If the residue is at either end of the protein sequence, an 'end move' is applied.
        - If the residue is elsewhere, one of 'corner move', 'crankshaft move', or 'pull move' is randomly chosen based on provided weights and applied.
        The probability weights of the moves depend on the `pm_weight`. The 'corner move' and 'crankshaft move' each have a weight of (0.5 - pm_weight/2), whereas the 'pull move' has a weight of pm_weight.
        If the selected movement is unsuccessful (indicated by the boolean 'success'), a new residue index is randomly chosen, and the process is repeated.

        Parameters:
            residue_index (int): The index of the residue to which movement is to be applied.

        Returns:
            None
        """
        success = False
        while not success:
            # Check if the residue_index is at either end of the sequence
            if residue_index in [0, len(self.sequence) - 1]:
                # Apply end move for terminal residues
                success = self.end_move(residue_index)
            else:
                # Choose one of the three moves with given probabilities
                move_type = self.weighted_random_choice(['corner', 'crankshaft', 'pull'], [0.5 - pm_weight/2, 0.5 - pm_weight/2, pm_weight])  # all the weights sum up to 1

                # Apply the selected move
                if move_type == 'corner':
                    success = self.corner_move(residue_index)
                elif move_type == 'crankshaft':
                    success = self.crankshaft_move(residue_index)
                else:  # pull move
                    success = self.pull_move(residue_index)

            # If move failed, pick a new residue_index for the next loop iteration until a move is made
            if not success:
                residue_index = random.randint(0, len(self.sequence) - 1)

    def compute_energy(self):
        """Compute the hydrophobic interaction energy of a conformation.
        This method quantifies the energy based on Dill's HP model, considering the interactions between hydrophobic (H) residues. 
        An energy penalty of -1 is applied for each pair of hydrophobic residues that are topological neighbours, 
        i.e. hydrophobic residues that are adjacent spatially but not sequentially.

        Returns:
            float: The energy of the conformation.
        """
        energy = 0
        sequence_length = len(self.sequence)

        # Set to keep track of checked residue pairs to avoid double counting the hydrophobic interactions
        passed_neighbours = set()  # e.g., if residue 1 is adjacent to residue 2, then residue 2 is also adjacent to residue 1

        for residue_index in range(sequence_length):
            # Get the polarity classification of the residue at the given index
            residue = self.residues[residue_index].polarity_classification

            # If the residue is hydrophobic, check its topological neighbors
            if residue == "H":
                x, y = self.positions[residue_index]
                adj_positions = [[x-1, y], [x+1, y], [x, y-1], [x, y+1]]

                for position in adj_positions:
                    # Ensure the adjacent position is occupied by a residue
                    if position in self.positions:
                        neighbour_index = self.positions.index(position) 
                        neighbour_residue = self.residues[neighbour_index].polarity_classification 

                        # Check the conditions for counting the energy interaction:
                        # 1. The residues should not be sequential neighbors
                        # 2. The pair should not have been evaluated before
                        # 3. The neighbor should also be hydrophobic
                        conditions_met = abs(neighbour_index - residue_index) != 1
                        conditions_met &= (residue_index, neighbour_index) not in passed_neighbours
                        conditions_met &= (neighbour_index, residue_index) not in passed_neighbours
                        conditions_met &= neighbour_residue == "H"

                        # If all conditions are met, update the energy and mark the pair as evaluated
                        if conditions_met:
                            energy -= 1
                            passed_neighbours.add((residue_index, neighbour_index))
                            passed_neighbours.add((neighbour_index, residue_index))

        # Store the computed energy in the class attribute
        self.energy = energy
        return self.energy
