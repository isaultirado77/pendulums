import numpy as np
from simulation import Mass

class SimplePendulum:
    def __init__(self, mass, length, theta0, omega0, g=9.81):
        """
        Inicializa un péndulo simple.
        :param mass: Masa del péndulo.
        :param length: Longitud del brazo.
        :param theta0: Ángulo inicial en radianes.
        :param omega0: Velocidad angular inicial en rad/s.
        :param g: Aceleración gravitacional (m/s^2).
        """
        self.mass = Mass(mass, length, theta0, omega0)
        self.g = g

    def equations_of_motion(self, state):
        """Devuelve las ecuaciones de movimiento del péndulo simple."""
        theta, omega = state[:2]
        dtheta_dt = omega
        domega_dt = -(self.g / self.mass.length) * np.sin(theta)
        return np.array([dtheta_dt, domega_dt])

    def E_kin(self):
        """Calcula la energía cinética del péndulo."""
        return 0.5 * self.mass.mass * (self.mass.length * self.mass.omega) ** 2

    def E_pot(self):
        """Calcula la energía potencial del péndulo."""
        return self.mass.mass * self.g * self.mass.length * (1 - np.cos(self.mass.theta))

    def E_tot(self):
        """Calcula la energía total del péndulo."""
        return self.E_kin() + self.E_pot()

