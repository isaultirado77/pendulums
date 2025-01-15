import numpy as np

class Pendulum:
    def __init__(self, length, mass, damping, driving_amplitude, driving_frequency):
        self.length = length
        self.mass = mass
        self.damping = damping
        self.driving_amplitude = driving_amplitude
        self.driving_frequency = driving_frequency
        self.theta = 0.2  # Initial angle (rad)
        self.omega = 0.0  # Initial angular velocity (rad/s)
        self.g = 9.81  # Gravitational acceleration (m/s^2)

    def f1(self, t, theta, omega):
        """First derivative of theta (angular velocity)."""
        return omega

    def f2(self, t, theta, omega):
        """Second derivative of theta (angular acceleration)."""
        term1 = self.driving_amplitude / (self.length * self.mass) * np.sin(self.driving_frequency * t)
        term2 = -(self.g / self.length) * np.sin(theta)
        term3 = -(self.damping / self.mass) * omega
        return term1 + term2 + term3

    def rk4_step(self, t, dt):
        """
        Perform a single step of the RK4 method.

        :param t: Current time (s)
        :param dt: Time step (s)
        """
        k1_theta = self.f1(t, self.theta, self.omega)
        k1_omega = self.f2(t, self.theta, self.omega)

        k2_theta = self.f1(t + 0.5 * dt, self.theta + 0.5 * k1_theta * dt, self.omega + 0.5 * k1_omega * dt)
        k2_omega = self.f2(t + 0.5 * dt, self.theta + 0.5 * k1_theta * dt, self.omega + 0.5 * k1_omega * dt)

        k3_theta = self.f1(t + 0.5 * dt, self.theta + 0.5 * k2_theta * dt, self.omega + 0.5 * k2_omega * dt)
        k3_omega = self.f2(t + 0.5 * dt, self.theta + 0.5 * k2_theta * dt, self.omega + 0.5 * k2_omega * dt)

        k4_theta = self.f1(t + dt, self.theta + k3_theta * dt, self.omega + k3_omega * dt)
        k4_omega = self.f2(t + dt, self.theta + k3_theta * dt, self.omega + k3_omega * dt)

        self.theta += (dt / 6) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta)
        self.omega += (dt / 6) * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega)

    def simulate(self, total_time, time_step):
        n_steps = int(total_time / time_step)
        time = np.linspace(0, total_time, n_steps)
        theta_values = []
        omega_values = []
        kinetic_energy = []
        potential_energy = []
        total_energy = []

        for t in time:
            theta_values.append(self.theta)
            omega_values.append(self.omega)

            ke = 0.5 * self.mass * (self.length * self.omega) ** 2
            pe = self.mass * self.g * self.length * (1 - np.cos(self.theta))
            kinetic_energy.append(ke)
            potential_energy.append(pe)
            total_energy.append(ke + pe)

            self.rk4_step(t, time_step)

        return time, np.array(theta_values), np.array(omega_values), np.array(kinetic_energy), np.array(potential_energy), np.array(total_energy)

# Example usage
if __name__ == "__main__":
    pendulum = Pendulum(length=9.8, mass=1.0, damping=0.5, driving_amplitude=0.5, driving_frequency=23.0)
    total_time = 150  # seconds
    time_step = 0.04  # seconds

    time, theta, omega, ke, pe, te = pendulum.simulate(total_time, time_step)

    # Save results to a file
    np.savetxt("pendulum_simulation.csv", np.column_stack([time, theta, omega, ke, pe, te]),
               header="Time,Theta,Omega,Kinetic Energy,Potential Energy,Total Energy", delimiter=",", comments="")
