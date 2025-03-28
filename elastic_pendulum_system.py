import os, argparse
import numpy as np

class ElasticPendulum:
    def __init__(self, mass, l0, k, theta0, omega0, l_init, v_init, g=9.81):
        """Inicializa un péndulo elástico."""
        self.mass = mass
        self.l0 = l0  # Longitud natural del resorte
        self.k = k    # Constante elástica
        self.g = g
        self.theta = theta0
        self.omega = omega0
        self.l = l_init
        self.v = v_init  # Velocidad inicial de la longitud del resorte

    def equations_of_motion(self, state):
        """Ecuaciones de movimiento del péndulo elástico."""
        theta, omega, l, v = state
        dtheta_dt = omega
        domega_dt = - (self.g / l) * np.sin(theta) - 2 * (v / l) * omega
        dl_dt = v
        dv_dt = - (self.k / self.mass) * (l - self.l0) + self.g * np.cos(theta) + l * omega ** 2
        return np.array([dtheta_dt, domega_dt, dl_dt, dv_dt])

def runge_kutta_step(pendulum, state, dt):
    """Realiza un paso de integración usando el método de Runge-Kutta de 4to orden."""
    k1 = dt * pendulum.equations_of_motion(state)
    k2 = dt * pendulum.equations_of_motion(state + 0.5 * k1)
    k3 = dt * pendulum.equations_of_motion(state + 0.5 * k2)
    k4 = dt * pendulum.equations_of_motion(state + k3)
    return state + (k1 + 2*k2 + 2*k3 + k4) / 6

def simulate_elastic_pendulum(pendulum, t_max, dt, filename="elastic_pendulum_data"):
    """Ejecuta la simulación del péndulo elástico y guarda los datos en un archivo."""
    os.makedirs("data", exist_ok=True)
    filepath = f"data/{filename}.dat"

    steps = int(t_max / dt)
    state = np.array([pendulum.theta, pendulum.omega, pendulum.l, pendulum.v])

    with open(filepath, "w") as file:
        file.write("# t theta omega l v E_kin E_pot E_elastic E_tot\n")
        time = 0.0
        for _ in range(steps):
            theta, omega, l, v = state
            E_kin = 0.5 * pendulum.mass * ((l * omega) ** 2 + v ** 2)
            E_pot = pendulum.mass * pendulum.g * l * (1 - np.cos(theta))
            E_elastic = 0.5 * pendulum.k * (l - pendulum.l0) ** 2
            E_tot = E_kin + E_pot + E_elastic

            file.write(f"{time:.5f} {theta:.5f} {omega:.5f} {l:.5f} {v:.5f} {E_kin:.5f} {E_pot:.5f} {E_elastic:.5f} {E_tot:.5f}\n")

            state = runge_kutta_step(pendulum, state, dt)
            time += dt

    print(f"Simulación completada. Datos guardados en {filepath}")

def main(mass=1.0, l0=1.0, k=10.0, theta0=np.pi/4, omega0=0.1, l_init=1.2, v_init=0.0, dt=0.01, t_max=10.0, filename="elastic_pendulum"):
    pendulum = ElasticPendulum(mass, l0, k, theta0, omega0, l_init, v_init)
    simulate_elastic_pendulum(pendulum, t_max, dt, filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulador de péndulo elástico usando Runge-Kutta de 4to orden.")
    parser.add_argument("--mass", type=float, default=1.0, help="Masa del péndulo (kg)")
    parser.add_argument("--l0", type=float, default=1.0, help="Longitud natural del resorte (m)")
    parser.add_argument("--k", type=float, default=10.0, help="Constante elástica del resorte (N/m)")
    parser.add_argument("--theta0", type=float, default=np.pi/4, help="Ángulo inicial (radianes)")
    parser.add_argument("--omega0", type=float, default=0.1, help="Velocidad angular inicial (rad/s)")
    parser.add_argument("--l_init", type=float, default=1.2, help="Longitud inicial del resorte (m)")
    parser.add_argument("--v_init", type=float, default=0.0, help="Velocidad inicial de la longitud del resorte (m/s)")
    parser.add_argument("--dt", type=float, default=0.01, help="Paso de tiempo para la simulación (s)")
    parser.add_argument("--t_max", type=float, default=10.0, help="Tiempo total de simulación (s)")
    parser.add_argument("--filename", type=str, default="elastic_pendulum", help="Nombre del archivo de salida")
    
    args = parser.parse_args()
    main(args.mass, args.l0, args.k, args.theta0, args.omega0, args.l_init, args.v_init, args.dt, args.t_max, args.filename)
