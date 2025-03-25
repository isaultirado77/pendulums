import numpy as np
from simutils import Mass


class ForcedPendulum:
    def __init__(self, mass, length, theta0, omega0, damping=0.1, driving_force=1.0, driving_freq=1.0, g=9.81):
        """
        Inicializa un péndulo forzado.
        :param mass: Masa del péndulo.
        :param length: Longitud del brazo.
        :param theta0: Ángulo inicial en radianes.
        :param omega0: Velocidad angular inicial en rad/s.
        :param damping: Coeficiente de amortiguamiento.
        :param driving_force: Amplitud de la fuerza externa.
        :param driving_freq: Frecuencia de la fuerza externa.
        :param g: Aceleración gravitacional (m/s^2).
        """
        self.mass = Mass(mass, length, theta0, omega0)
        self.g = g
        self.damping = damping
        self.driving_force = driving_force
        self.driving_freq = driving_freq

    def equations_of_motion(self, state, time):
        """Devuelve las ecuaciones de movimiento del péndulo forzado."""
        theta, omega = state
        dtheta_dt = omega
        domega_dt = -(self.g / self.mass.length) * np.sin(theta) - self.damping * omega + self.driving_force * np.cos(self.driving_freq * time)
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
