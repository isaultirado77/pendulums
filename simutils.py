class Mass:
    def __init__(self, mass, length, theta, omega):
        """
        Inicializa una masa de un péndulo.
        :param mass: Masa del péndulo.
        :param length: Longitud del brazo.
        :param theta: Ángulo inicial (en radianes).
        :param omega: Velocidad angular inicial (en rad/s).
        """
        self.mass = mass
        self.length = length
        self.theta = theta
        self.omega = omega

