import pandas as pd
from beamline import *


class ExcelElements:
    def __init__(self, file_path):
        """
        Initialize the ExcelElements class and load the Excel file.

        :param file_path: Path to the Excel file containing the beamline information.
        """
        self.df = pd.DataFrame()  # DataFrame to store the entire spreadsheet data
        self.load_excel_lattice(file_path)

    def load_excel_lattice(self, file_path):
        """
        Load the lattice from an Excel file and store it in a DataFrame.

        :param file_path: Path to the Excel file.
        """
        # Load the Excel file without headers and skip the first row
        df = pd.read_excel(file_path, header=None, skiprows=1)

        # Rename columns to match their usage
        df.columns = ['Nomenclature', 'z_sta', 'z_mid', 'z_end', 'I (A)', 'Element name', 'ch', 'Difference', 'Sector',...,
                      'Element','Vacuum sectors','z logbook UH (m)']

        # Convert 'ch' to numeric, handling any non-numeric gracefully
        df['ch'] = pd.to_numeric(df['ch'], errors='coerce')

        # Store the DataFrame for external access
        self.df = df

    def create_beamline(self):
        """
        Create the beamline by iterating through the DataFrame and generating
        the corresponding beamline elements (QPF, QPD, drifts).

        :return: List of beamline elements (drifts, quadrupoles).
        """
        beamline = []  # This will hold the beamline elements
        prev_z_end = 0  # Track the end position of the previous element
        default_current = 0

        for index, row in self.df.iterrows():
            element = row['Element']
            z_sta = row['z_sta']
            z_end = row['z_end']

            # Read the current for the quadrupole; if not a number, default to 0
            try:
                current = float(row['I (A)']) if pd.notna(row['I (A)']) else 0.0
            except (ValueError, TypeError):
                current = default_current

            # Calculate the drift length between previous element and current element
            if z_sta > prev_z_end:
                drift_length = z_sta - prev_z_end
                beamline.append(driftLattice(drift_length))

            # Add quadrupole focusing (QPF) or defocusing (QPD) based on nomenclature
            if element == "QPF":
                beamline.append(qpfLattice(current=current))
            elif element == "QPD":
                beamline.append(qpdLattice(current=current))
            else:
                # Add additional elements here if necessary
                pass

            # Update prev_z_end for the next element
            prev_z_end = z_end

        return beamline

    def get_dataframe(self):
        """
        Return the DataFrame containing the beamline elements.

        :return: DataFrame with the loaded Excel data.
        """
        return self.df

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
