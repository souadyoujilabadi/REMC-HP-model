"""A module that represents protein residues.

This module provides the `Residue` class which models a protein residue.
The class offers functionality to determine the polarity of a given protein
residue using Dill's Hydrophobic-Polar (HP) model.

Classes:
    - Residue: Represents a protein residue and provides a method to determine its polarity.

Example:
    >>> res = Residue("V")
    >>> print(res.polarity_classification)
    H
"""


class Residue:
    """Represents a protein residue.

    This class models a protein residue and provides a conversion mechanism according Dill's HP (Hydrophobic Polar) model.

    Attributes:
        residue (str): The residue letter.
        polarity_classification (str): The residue's polarity classification (H or P).

    Methods:
        convert_to_hp(): Convert the residue based on Dill's HP model.

    """

    POLARITY_MAP = {
        "V": "H", "I": "H", "F": "H", "L": "H", "M": "H", "C": "H", "W": "H",
        "A": "P", "Y": "P", "G": "P", "P": "P", "T": "P", "K": "P", "R": "P",
        "H": "P", "D": "P", "E": "P", "S": "P", "N": "P", "Q": "P", "U": "P"
    }

    def __init__(self, residue):
        """Initialize a Residue object with a given residue letter.

        Parameters:
            residue (str): A single-character string representing the residue.
        """
        self.residue = residue
        self.polarity_classification = self.convert_to_hp()

    def convert_to_hp(self):
        """Convert the residue based on Dill's HP (Hydrophobic Polar) model.

        Returns:
            str: The residue's polarity classification (H for Hydrophobic, P for Polar).
        """
        # Use the get() method of the dictionary to retrieve the polarity classification
        # for the current residue. If the residue is not found in the dictionary,
        # it defaults to returning "Unknown".
        return self.POLARITY_MAP.get(self.residue, "Unknown")