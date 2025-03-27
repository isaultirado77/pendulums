import os
import numpy as np

class ForcedPendulum:
    def __init__(self, mass, length, theta0, omega0, A, omega_f, g=9.81):
        """Inicializa un péndulo forzado."""
        self.mass = mass
        self.length = length
        self.theta = theta0
        self.omega = omega0
        self.A = A 
        self.omega_f = omega_f 
        self.g = g

    def equations_of_motion(self, state, t):
        """Ecuaciones de movimiento del péndulo forzado."""
        theta, omega = state
        dtheta_dt = omega
        domega_dt = -(self.g / self.length) * np.sin(theta) + self.A * np.cos(self.omega_f * t)
        return np.array([dtheta_dt, domega_dt])

def runge_kutta_step(pendulum, state, t, dt):
    """Realiza un paso de integración usando el método de Runge-Kutta de 4to orden."""
    k1 = dt * pendulum.equations_of_motion(state, t)
    k2 = dt * pendulum.equations_of_motion(state + 0.5 * k1, t + 0.5 * dt)
    k3 = dt * pendulum.equations_of_motion(state + 0.5 * k2, t + 0.5 * dt)
    k4 = dt * pendulum.equations_of_motion(state + k3, t + dt)
    return state + (k1 + 2*k2 + 2*k3 + k4) / 6

def simulate_forced_pendulum(pendulum, t_max, dt, filename="forced_pendulum_data"):
    """Ejecuta la simulación del péndulo forzado y guarda los datos en un archivo."""
    os.makedirs("data", exist_ok=True)
    filepath = f"data/{filename}.dat"
    
    steps = int(t_max / dt)
    state = np.array([pendulum.theta, pendulum.omega])
    
    with open(filepath, "w") as file:
        file.write("# t theta omega E_kin E_pot E_tot\n")
        time = 0.0
        for _ in range(steps):
            theta, omega = state
            E_kin = 0.5 * pendulum.mass * (pendulum.length * omega) ** 2
            E_pot = pendulum.mass * pendulum.g * pendulum.length * (1 - np.cos(theta))
            E_tot = E_kin + E_pot

            file.write(f"{time:.5f} {theta:.5f} {omega:.5f} {E_kin:.5f} {E_pot:.5f} {E_tot:.5f}\n")

            state = runge_kutta_step(pendulum, state, time, dt)
            time += dt

    print(f"Simulación completada. Datos guardados en {filepath}")

def main():
    # Parámetros del sistema
    mass = 1.0
    length = 1.0
    theta0 = np.pi / 4
    omega0 = 0.1

    # Parámetros de la fuerza forzada
    A = 1.0   # Amplitud de la fuerza externa
    omega_f = 1.5  # Frecuencia de la fuerza externa

    # Parámetros de simulación
    dt = 0.01
    t_max = 10.0
    filename = "test_forced_pendulum"

    pendulum = ForcedPendulum(mass, length, theta0, omega0, A, omega_f)
    simulate_forced_pendulum(pendulum, t_max, dt, filename)

if __name__ == "__main__":
    main()
