import pandas as pd
import numpy as np

#should I use the position that each beamline object takes up instead of its total length?

class lattice:
    def __init__(self, length: float = -1):
        self.nomenclatures = []  # List of strings from the first column
        self.positions = pd.DataFrame()  # DataFrame for numerical columns
        self.descriptions = []  # List of descriptions from the sixth column
        self.ch = []  # List (array) of channel numbers from the last column

        self.E = 35  # Kinetic energy (MeV/c^2)
        self.E0 = 0.51099
        self.length = length
        self.gamma = (1 + (self.E/self.E0))

    
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
            # Check if z is within the range z_sta to z_end
            if row['z_sta'] <= z <= row['z_end']:
                # Retrieve the corresponding nomenclature
                nomenclature = self.nomenclatures[index]
                # You can customize this part to perform specific string comparisons or outputs
                # For simplicity, let's return the nomenclature
                return nomenclature
        # If z does not fall within any range, return a default value
        return "Position out of range"

    def __str__(self):
        return f"Beamline: {len(self.nomenclatures)} elements"
    



class driftLattice(lattice):
    def __init__(self, length: float = -1):
        if length != -1:
            super().__init__(length = length)
        else: raise ValueError("Invalid Parameter: Please enter a positive length parameter")

    #Matrix multiplecation, values is a 2 dimensional numPy array, each array is 6 elements long
    #values = np.array([[x, x', y, y', z, z'],...])
    #Note: 1x6 array is multiplied correctly with 6x6 array
    def getDriftMatrice(self, values, length = -1):
        if length == -1:
            length = self.length

        matrice = (np.array([[1, length, 0, 0, 0, 0],
                                 [0, 1, 0, 0, 0, 0],
                                 [0, 0, 1, length, 0, 0],
                                 [0, 0, 0, 1, 0, 0],
                                 [0, 0, 0, 0, 1, (length/((self.gamma)**2))],
                                 [0, 0, 0, 0, 0, 1]]))
        newMatrix = []
        for array in values:
            tempArray = np.matmul(matrice, array)
            newMatrix.append(tempArray.tolist())
        return newMatrix #  return 2d list






class qpfLattice(lattice):
    def __init__(self, length: float = -1, current: float = 0):
        self.current = current
        self.theta = np.sqrt(self.current)*length
        if length == -1:
            super().__init__(length = 88.9)
        else:
            super().__init__(length = length)

    '''
    performs a transformation to a 2d np array made of 1x6 variable matrices

    values: np.array([list[int],...])
    '''
    def getQPFmatrice(self, values, length = -1):
        if length == -1:
            length = self.length
        
        matrice = (np.array([[np.cos(self.theta),(np.sin(self.theta)/np.sqrt(self.current)),0,0,0,0],
                               [(-(np.sqrt(self.current)))*(np.sin(self.theta)),np.cos(self.theta),0,0,0,0],
                               [0,0,np.cosh(self.theta),(np.sinh(self.theta))/(np.sqrt(self.current)),0,0],
                               [0,0,np.sqrt(self.current)*np.sinh(self.theta),np.cosh(self.theta),0,0],
                               [0,0,0,0,1,length/((self.gamma)**2)],
                               [0,0,0,0,0,1]]))
        
        newMatrix = []
        for array in values:
            tempArray = np.matmul(matrice, array)
            newMatrix.append(tempArray.tolist())
        return newMatrix