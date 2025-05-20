import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

class Radiation:
    E_e_MeV = 45                    # Electron energy [MeV]
    e = 1.602176634e-19  # Elementary charge (C) [CODATA 2018]
    MeV_to_J = 1.602176634e-13  # Conversion factor from MeV to Joules [CODATA 2018]
    NA = 6.02214076e23  # Avogadro's number
    epsilon_0 = 8.854187817e-12  # Vacuum permittivity (F/m)
    me = 9.10938356e-31  # Electron mass (kg)
    m_p = 1.67262192595e-27  # Proton Mass (kg)
    c = 299792458.0  # Speed of light (m/s)
    r_e = e**2 / (4 * np.pi * epsilon_0 * me * c**2)  # Classical electron radius [m]
    h = 6.62607015e-34        # Planck constant (J·s)
    me_c2_J = me * c**2             # Electron rest energy [J]
    gamma = E_e_MeV * 1e6 * e  / me_c2_J         # Lorentz factor

    def __init__(self) -> None:
        pass

    # E_thresh_eV : 10 keV
    # lambda_L_um : Laser wavelength [um]
    def plot_ICS_angularDist(self, lambda_L_um = 3,E_thresh_eV = 1e4):

        # Threshold for desired photon energy [eV]
        E_thresh_J = E_thresh_eV * self.e

        lambda_L = lambda_L_um * 1e-6                  # Laser wavelength [m]
        E_gamma_L = self.h * self.c / lambda_L                   # Photon energy [J]

        # Find theta_max where E_gamma(theta) = E_thresh
        # For relativistic beams, approximate small angle solution:
        theta_thresh_rad = np.sqrt((4 * self.gamma**2 * E_gamma_L / E_thresh_J - 1 - 4 * self.gamma * E_gamma_L / self.me_c2_J) / self.gamma**2)
        theta_thresh_mrad = theta_thresh_rad * 1e3

        
        
        theta_vals = np.linspace(0, 5/self.gamma, 1000)
        theta_mrad = theta_vals * 1e3

        
        sigma_T = (8 * np.pi / 3) * self.r_e**2  # # Thomson scattering cross section [m^2]
        # Compute differential cross section (Thomson) * sin(theta)
        d_sigma_vals = 0.5 * self.r_e**2 * (1 + np.cos(theta_vals)**2) * np.sin(theta_vals)
        # Normalize to total cross section
        d_sigma_normalized = d_sigma_vals / sigma_T

        # Plot
        plt.figure(figsize=(8,5))
        plt.plot(theta_mrad, d_sigma_normalized, label=r"$\frac{1}{\sigma_T}\frac{d\sigma}{d\Omega} \cdot \sin\theta$", color='blue')
        plt.axvline(theta_thresh_mrad, color='red', linestyle='--', label="10 keV cutoff")
        plt.fill_between(theta_mrad, d_sigma_normalized, where=(theta_mrad <= theta_thresh_mrad), color='blue', alpha=0.3, label="E ≥ 10 keV")
        plt.xlabel("Scattering angle θ [mrad]")
        plt.ylabel("Normalized angular distribution")
        plt.title("ICS Angular Distribution (to be fixed)")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    obj = Radiation()
    obj.plot_ICS_angularDist()