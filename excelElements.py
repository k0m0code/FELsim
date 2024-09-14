import pandas as pd
import numpy as np
from beamline import *

class excelElements:
    def get_matrix(self, identifier):
        # Split the identifier by dots and extract the middle part
        parts = identifier.split('.')
        if len(parts) != 3:
            raise ValueError("Identifier must be in the format 'XXX.YYY.123'")

        middle_part = parts[1]

        # Define matrices for each case
        matrices = {
            "QPF": np.array([[1, 0, 0, 0, 0, 0],  # Placeholder matrix for QPF
                            [0, 1, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0]
                            ]),
            "QPD": np.array([[2, 0, 0, 0, 0, 0],  # Placeholder matrix for QPD
                            [0, 2, 0, 0, 0, 0],
                            # Add the rest of the rows for a 6x6 matrix
                            ]),
            # Add matrices for SB1, SB2, ... SB5 in a similar manner
            "SB1": np.array([[1, 0, 0, 0, 0, 0],  # Placeholder matrix for SB1
                            [0, 1, 0, 0, 0, 0],
                            # Add the rest of the rows for a 6x6 matrix
                            ]),
            # Continue for SB2 to SB5...
        }

        # Check if the middle part matches any key in the matrices dictionary
        if middle_part in matrices:
            return matrices[middle_part]
        else:
            raise ValueError(f"No matrix defined for identifier {middle_part}")

    def load_excel_lattice(self, file_path):
        # Load the Excel file without headers and skip the first row
        df = pd.read_excel(file_path, header=None, skiprows=1)

        # Assuming the first column contains string identifiers
        self.nomenclatures = df.iloc[:, 0].tolist()

        # Columns 2 to 5 contain numbers, creating a DataFrame directly
        self.positions = df.iloc[:, 1:5]
        self.positions.columns = ['z_sta', 'z_mid', 'z_end', 'z_log']

        # Column 6 is a description
        self.descriptions = df.iloc[:, 5].tolist()

        # The last column contains channel numbers
        self.ch = df.iloc[:, -1].tolist()

        # Convert 'ch' to a numeric type, handling any non-numeric gracefully
        self.ch = pd.to_numeric(self.ch, errors='coerce').tolist()

    def find_element_by_position(self, z):
        # Iterate through the DataFrame rows
        for index, row in self.positions.iterrows():
            # Check if z is within the range  z_sta to z_end
            if row['z_sta'] <= z <= row['z_end']:
                # Retrieve the corresponding nomenclature
                nomenclature = self.nomenclatures[index]
                # You can customize this part to perform specific string comparisons or outputs
                # For simplicity, let's return the nomenclature
                return nomenclature
        # If z does not fall within any range, return a default value
        return "Position out of range"
    
    def loadBeamSegments(self):
        beamline = []
        for index, row in self.positions.iterrows():
            nomenclature = self.nomenclatures[index]
            zend = row['z_end']
            zstart = row['z_sta']
            length = zend-zstart
            if nomenclature.upper() == "QPF":
                segment = qpfLattice()


    def __str__(self):
        return f"Beamline: {len(self.nomenclatures)} elements"
