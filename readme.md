# Simulación de Péndulo Físico

Este proyecto simula el movimiento de un péndulo físico bajo la influencia de la gravedad, amortiguamiento y una fuerza impulsora, utilizando Python. La simulación emplea el método de Runge-Kutta de cuarto orden (RK4) para resolver las ecuaciones de movimiento.

## Características

- Simula un péndulo forzado con amortiguamiento.
- Permite personalizar parámetros del péndulo como la longitud, masa, coeficiente de amortiguamiento, amplitud de la fuerza impulsora y frecuencia angular de la fuerza impulsora.
- Calcula y almacena el desplazamiento angular (theta), velocidad angular (omega), energía cinética, energía potencial y energía total a lo largo del tiempo.
- Guarda los resultados de la simulación en un archivo CSV para análisis posterior.

## Requisitos Previos

- Python 3.8 o superior.
- Librerías necesarias:
  - `numpy`
  - `matplotlib`

Instala las librerías requeridas usando pip:
```bash
pip install numpy matplotlib
```

## Ejecución de la Simulación

1. Clona el repositorio o descarga el script.
2. Ejecuta el script:
   ```bash
   python pendulum_simulation.py
   ```
3. Los resultados de la simulación se guardarán en el archivo `pendulum_simulation.csv` en el directorio data.

### Ejemplo de Uso

El script incluye un ejemplo de uso de la clase `Pendulum`. Simula un péndulo con los siguientes parámetros:

- Longitud: 9.8 m
- Masa: 1.0 kg
- Coeficiente de amortiguamiento: 0.5 kg/s
- Amplitud de la fuerza impulsora: 0.5 N
- Frecuencia angular de la fuerza impulsora: 23 rad/s
- Tiempo de simulación: 150 segundos
- Paso de tiempo: 0.04 segundos

## Salida

La simulación genera un archivo CSV con las siguientes columnas:

- `Time`: Tiempo (s).
- `Theta`: Desplazamiento angular (rad).
- `Omega`: Velocidad angular (rad/s).
- `Kinetic Energy`: Energía cinética (J).
- `Potential Energy`: Energía potencial (J).
- `Total Energy`: Energía total (J).

## Visualización

Puedes visualizar los resultados utilizando cualquier herramienta de gráficos, como Matplotlib. Por ejemplo:

```python
import matplotlib.pyplot as plt
import numpy as np

# Cargar los datos
data = np.loadtxt("pendulum_simulation.dat", delimiter=",", skiprows=1)

# Extraer columnas
tiempo = data[:, 0]
theta = data[:, 1]
omega = data[:, 2]
energia_cinetica = data[:, 3]
energia_potencial = data[:, 4]
energia_total = data[:, 5]

# Graficar desplazamiento angular
plt.figure()
plt.plot(tiempo, theta, label="Theta (rad)")
plt.xlabel("Tiempo (s)")
plt.ylabel("Desplazamiento Angular (rad)")
plt.title("Desplazamiento Angular vs Tiempo")
plt.legend()
plt.show()
```
---