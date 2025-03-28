import os
import numpy as np
import argparse

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
        self.mass = mass
        self.length = length
        self.theta = theta0
        self.omega = omega0
        self.g = g

    def equations_of_motion(self, state):
        """Devuelve las ecuaciones de movimiento del péndulo simple."""
        theta, omega = state
        dtheta_dt = omega
        domega_dt = -(self.g / self.length) * np.sin(theta)
        return np.array([dtheta_dt, domega_dt])

    def E_kin(self):
        """Calcula la energía cinética del péndulo."""
        return 0.5 * self.mass * (self.length * self.omega) ** 2

    def E_pot(self):
        """Calcula la energía potencial del péndulo."""
        return self.mass * self.g * self.length * (1 - np.cos(self.theta))

    def E_tot(self):
        """Calcula la energía total del péndulo."""
        return self.E_kin() + self.E_pot()

class Simulator:
    def __init__(self, pendulum, filename="pendulum_data"):
        """
        Inicializa el sistema de simulación del péndulo.
        :param pendulum: Instancia de un péndulo específico.
        :param filename: Nombre del archivo donde se guardarán los datos.
        """
        self.pendulum = pendulum
        self.filename = filename
        os.makedirs("data", exist_ok=True)

    def runge_kutta_step(self, state, dt):
        """Realiza un paso de integración usando el método de Runge-Kutta de 4to orden."""
        k1 = dt * self.pendulum.equations_of_motion(state)
        k2 = dt * self.pendulum.equations_of_motion(state + 0.5 * k1)
        k3 = dt * self.pendulum.equations_of_motion(state + 0.5 * k2)
        k4 = dt * self.pendulum.equations_of_motion(state + k3)
        return state + (k1 + 2*k2 + 2*k3 + k4) / 6

    def simulate(self, t_max, dt):
        """Ejecuta la simulación y guarda los datos en un archivo."""
        steps = int(t_max / dt)
        state = np.array([self.pendulum.theta, self.pendulum.omega])
        
        with open(f"data/{self.filename}.dat", "w") as file:
            file.write("# t theta omega E_kin E_pot E_tot\n")
            time = 0.0
            for _ in range(steps):
                theta, omega = state
                self.pendulum.theta = theta  # Actualizar valores actuales
                self.pendulum.omega = omega
                E_kin = self.pendulum.E_kin()
                E_pot = self.pendulum.E_pot()
                E_tot = self.pendulum.E_tot()
                file.write(f"{time:.5f} {theta:.5f} {omega:.5f} {E_kin:.5f} {E_pot:.5f} {E_tot:.5f}\n")
                state = self.runge_kutta_step(state, dt)
                time += dt

        print(f"Simulación completada. Datos guardados en data/{self.filename}.dat")


def main(mass=1.0, length=1.0, theta0=np.pi/4, omega0=0.1, dt=0.01, t_max=10.0, filename="test_simple_pendulum"):
    pendulum = SimplePendulum(mass, length, theta0, omega0)
    simulator = Simulator(pendulum, filename)
    simulator.simulate(t_max, dt)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulador de péndulo simple usando Runge-Kutta de 4to orden.")
    parser.add_argument("--mass", type=float, default=1.0, help="Masa del péndulo (kg)")
    parser.add_argument("--length", type=float, default=1.0, help="Longitud del péndulo (m)")
    parser.add_argument("--theta0", type=float, default=np.pi/4, help="Ángulo inicial (radianes)")
    parser.add_argument("--omega0", type=float, default=0.1, help="Velocidad angular inicial (rad/s)")
    parser.add_argument("--dt", type=float, default=0.01, help="Paso de tiempo para la simulación (s)")
    parser.add_argument("--t_max", type=float, default=10.0, help="Tiempo total de simulación (s)")
    parser.add_argument("--filename", type=str, default="test_simple_pendulum", help="Nombre del archivo de salida")
    
    args = parser.parse_args()
    main(args.mass, args.length, args.theta0, args.omega0, args.dt, args.t_max, args.filename)
