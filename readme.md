
# Simulación de Péndulos

Este proyecto simula diferentes tipos de péndulos: simple, forzado-amortiguado y elástico, utilizando el método de Runge-Kutta de 4to orden para resolver las ecuaciones de movimiento. Se analizan trayectorias, espacio fase, evolución angular y energía.

## Análisis y Visualización

Los resultados de la simulación pueden explorarse en el siguiente notebook:

[![Ver Análisis en Jupyter](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.org/github/isaultirado77/pendulums/blob/main/pendulum_simulation_analysis.ipynb)


## Estructura del Proyecto

- **`data/`** – Almacena los datos generados por las simulaciones.
- **`plots/`** – Contiene gráficos de los resultados.
- **`animations/`** – Guarda animaciones en formato GIF.
- **`pendulum_simulation_analysis.ipynb`** – Notebook donde se analizan los datos generados.
- **Scripts principales**:
  - `simple_pendulum.py`
  - `forced_pendulum_system.py`
  - `elastic_pendulum_system.py`
  - `animator.py`

## Ejecución de la Simulación

Cada script puede ejecutarse con valores por defecto o configurarse con `argparse`.

Ejemplo para el péndulo forzado-amortiguado:

```bash
python forced_pendulum_system.py --mass 1.0 --length 1.0 --theta0 0.785 \
                                 --omega0 0.1 --A 1.0 --omega_f 1.5 --gamma 0.1 \
                                 --dt 0.01 --t_max 10.0 --filename output
```

## Requisitos

- **Python 3.7+**
- Instalación de dependencias:
  ```bash
  pip install numpy matplotlib scipy
  ```
