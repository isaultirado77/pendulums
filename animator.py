import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

class Animator:
    def __init__(self, t, theta, length=1.0, draw_trace=False):
        self.t = t
        self.theta = theta
        self.length = length
        self.draw_trace = draw_trace
        self.index = 0  # Índice para recorrer los datos

        # Convertir coordenadas angulares a cartesianas
        self.x, self.y = self.compute_cartesian_coords()

        # Configuración de la figura
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim(-1.5 * self.length, 1.5 * self.length)
        self.ax.set_ylim(-1.5 * self.length, 1.5 * self.length)

        # Texto para mostrar el tiempo
        self.time_text = self.ax.text(0.05, 0.95, '', 
            horizontalalignment='left', 
            verticalalignment='top', 
            transform=self.ax.transAxes)

        # Dibujar el punto de anclaje
        self.ax.plot(0, 0, 'ko', markersize=10)

        # Dibujar la masa y la línea del péndulo
        self.line, = self.ax.plot([], [], 'ro-', lw=2, markersize=8)  
        
        # Si draw_trace es True, dibuja la trayectoria de la masa
        if self.draw_trace:
            self.trace, = self.ax.plot([], [], 'r-', alpha=0.3)

    def compute_cartesian_coords(self):
        """Convierte coordenadas angulares a cartesianas."""
        x = self.length * np.sin(self.theta)
        y = -self.length * np.cos(self.theta)
        return x, y

    def advance_time_step(self):
        """Generador que recorre los datos precomputados."""
        while self.index < len(self.t):
            yield self.index
            self.index += 5

    def update(self, index):
        """Actualiza la animación en cada frame."""
        self.time_text.set_text(f'Tiempo: {self.t[index]:.2f} s')
        
        # Actualizar la posición del péndulo
        self.line.set_data([0, self.x[index]], [0, self.y[index]])

        # Actualizar la traza de la masa
        if self.draw_trace:
            self.trace.set_data(self.x[:index], self.y[:index])

        return self.line,

    def animate(self):
        """Ejecuta la animación."""
        self.animation = animation.FuncAnimation(
            self.fig, self.update, frames=self.advance_time_step,
            interval=5, blit=True
        )
    
    def save_animation(self, filename="single_pendulum.gif", fps=60):
        os.makedirs("animations", exist_ok=True)
        writer = animation.PillowWriter(fps=fps)
        self.animation.save(f"animations/{filename}", writer=writer)
