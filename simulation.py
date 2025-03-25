import numpy as np
import os
import argparse
import simple
import forced 
import damped
import elastic


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
        state = np.array([self.pendulum.mass.theta, self.pendulum.mass.omega])
        
        with open(f"data/{self.filename}.dat", "w") as file:
            file.write("# t theta omega E_kin E_pot E_tot\n")
            time = 0.0
            for _ in range(steps):
                theta, omega = state
                E_kin = self.pendulum.E_kin()
                E_pot = self.pendulum.E_pot()
                E_tot = self.pendulum.E_tot()
                file.write(f"{time:.5f} {theta:.5f} {omega:.5f} {E_kin:.5f} {E_pot:.5f} {E_tot:.5f}\n")
                state = self.runge_kutta_step(state, dt)
                time += dt

        print(f"Simulación completada. Datos guardados en data/{self.filename}.dat")

def main(): 
    # System parameters
    mass = 1
    length = 1

    # Simulation parameters
    dt = 0.01
    t_max = 10.0

    theta0 = np.pi / 4
    omega0 = 0.1
    simple_pendulum = simple.SimplePendulum(mass, length, theta0, omega0)
    filename = "simple_pendulum"

    sp_simulator = Simulator(simple_pendulum, filename)
    sp_simulator.simulate(t_max=t_max, dt=dt)

if __name__ == "__main__": 
    main()