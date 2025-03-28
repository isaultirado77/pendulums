import os
import numpy as np
import argparse

class ForcedDampedPendulum:
    def __init__(self, mass, length, theta0, omega0, A, omega_f, gamma, g=9.81):
        """Inicializa un péndulo forzado y amortiguado."""
        self.mass = mass
        self.length = length
        self.theta = theta0
        self.omega = omega0
        self.A = A  # Amplitud de la fuerza externa
        self.omega_f = omega_f  # Frecuencia de la fuerza externa
        self.gamma = gamma  # Coeficiente de amortiguamiento
        self.g = g

    def equations_of_motion(self, state, t):
        """Ecuaciones de movimiento del péndulo forzado y amortiguado."""
        theta, omega = state
        dtheta_dt = omega
        domega_dt = -(self.g / self.length) * np.sin(theta) - self.gamma * omega + self.A * np.cos(self.omega_f * t)
        return np.array([dtheta_dt, domega_dt])


def runge_kutta_step(pendulum, state, t, dt):
    """Realiza un paso de integración usando el método de Runge-Kutta de 4to orden."""
    k1 = dt * pendulum.equations_of_motion(state, t)
    k2 = dt * pendulum.equations_of_motion(state + 0.5 * k1, t + 0.5 * dt)
    k3 = dt * pendulum.equations_of_motion(state + 0.5 * k2, t + 0.5 * dt)
    k4 = dt * pendulum.equations_of_motion(state + k3, t + dt)
    return state + (k1 + 2*k2 + 2*k3 + k4) / 6


def simulate_forced_damped_pendulum(pendulum, t_max, dt, filename="forced_damped_pendulum_data"):
    """Ejecuta la simulación del péndulo forzado y amortiguado, guardando los datos en un archivo."""
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


def main(mass=1.0, length=1.0, theta0=np.pi/4, omega0=0.1, A=1.0, omega_f=1.5, gamma=0.1, dt=0.01, t_max=10.0, filename="test_forced_damped_pendulum"):
    pendulum = ForcedDampedPendulum(mass, length, theta0, omega0, A, omega_f, gamma)
    simulate_forced_damped_pendulum(pendulum, t_max, dt, filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulador de péndulo forzado y amortiguado usando Runge-Kutta de 4to orden.")
    parser.add_argument("--mass", type=float, default=1.0, help="Masa del péndulo (kg)")
    parser.add_argument("--length", type=float, default=1.0, help="Longitud del péndulo (m)")
    parser.add_argument("--theta0", type=float, default=np.pi/4, help="Ángulo inicial (radianes)")
    parser.add_argument("--omega0", type=float, default=0.1, help="Velocidad angular inicial (rad/s)")
    parser.add_argument("--A", type=float, default=1.0, help="Amplitud de la fuerza externa")
    parser.add_argument("--omega_f", type=float, default=1.5, help="Frecuencia de la fuerza externa (rad/s)")
    parser.add_argument("--gamma", type=float, default=0.1, help="Coeficiente de amortiguamiento")
    parser.add_argument("--dt", type=float, default=0.01, help="Paso de tiempo para la simulación (s)")
    parser.add_argument("--t_max", type=float, default=10.0, help="Tiempo total de simulación (s)")
    parser.add_argument("--filename", type=str, default="test_forced_damped_pendulum", help="Nombre del archivo de salida")
    
    args = parser.parse_args()
    main(args.mass, args.length, args.theta0, args.omega0, args.A, args.omega_f, args.gamma, args.dt, args.t_max, args.filename)
