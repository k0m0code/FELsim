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
        df.columns = ['Nomenclature', 'z_start', 'z_mid', 'z_end', 'Current (A)', 'Dipole Angle (deg)',
                      'Dipole length (m)', 'Dipole wedge (deg)', 'Gap wedge (m)', 'Element name',
                      'Channel','Label','Sector', 'Element']
        

        # Convert 'ch' to numeric, handling any non-numeric gracefully
        df['Channel'] = pd.to_numeric(df['Channel'], errors='coerce')

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
            z_sta = row['z_start']
            z_end = row['z_end']

            # Read the current for the quadrupole
            current = float(row['Current (A)']) if pd.notna(row['Current (A)']) else 0.0
            angle = float(row['Dipole Angle (deg)']) if pd.notna(row['Dipole Angle (deg)']) else 0.0
            curvature = float(row['Dipole length (m)']) if pd.notna(row['Dipole length (m)']) else 0.0
            angle_wedge = float(row['Dipole wedge (deg)']) if pd.notna(row['Dipole wedge (deg)']) else 0.0
            gap_wedge = float(row['Gap wedge (m)']) if pd.notna(row['Gap wedge (m)']) else 0.0

            # Calculate the drift length between previous element and current element
            if z_sta > prev_z_end:
                # print(str(prev_z_end) + " " + str(z_sta)) #testing
                drift_length = z_sta - prev_z_end
                beamline.append(driftLattice(drift_length))

            # Add quadrupole focusing (QPF) or defocusing (QPD) based on nomenclature
            # Possibly match the excel three letters for the element with the first three letters of the function to call
            # name param1, param2 if possible, dict optional params
            if element == "QPF":
                # print(str(z_sta) + " " + str(z_end)) #testing
                beamline.append(qpfLattice(current=current, length=(z_end - z_sta)))
            elif element == "QPD":
                # print(str(z_sta) + " " + str(z_end)) #testing
                beamline.append(qpdLattice(current=current, length=(z_end - z_sta)))
            elif element == "DPH":
                beamline.append(dipole(length=curvature, angle=angle))
            elif element == "DPW":
                beamline.append(dipole_wedge(length=gap_wedge, angle=angle_wedge, dipole_length=curvature,dipole_angle=angle))
            else:
                # # Add additional elements here if necessary
                # pass
                if (not z_end-z_sta == 0) and (not np.isnan(z_sta)) and (not np.isnan(z_end)):
                    beamline.append(driftLattice(z_end-z_sta))

            # Update prev_z_end for the next element
            if not np.isnan(z_end):
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


    def __str__(self):
        return f"Beamline: {len(self.nomenclatures)} elements"
