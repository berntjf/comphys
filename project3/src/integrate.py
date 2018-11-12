import numpy as np
from numba import jit

class IntegrationMethod:
    def __init__(self):
        self.xi = None
        self.vi = None
        self.ti = None

class ForwardEuler(IntegrationMethod):
    """
    Returnerer posisjon og hastighet sammen.
    """
    def __call__(self, f, x0, v0, t, steps=1):
        u = np.empty([t.size, 2] + list(x0.shape))
        u[0]      = np.array([x0, v0])
        current_u = u[0].copy()
        @jit
        def integrate(current_u):
            k = 0
            for k in range(1, t.size):
                dt = (t[k] - t[k-1]) / steps
                for i in range(steps):
                    current_u[0] = current_u[0] + dt*current_u[1]
                    current_u[1] = current_u[1] + dt*f(current_u[0])
                u[k] = current_u.copy()
            return u, k
        u, k = integrate(current_u)
        self.xi = u[:k+1,0,...]
        self.vi = u[:k+1,1,...]
        self.ti = t[:k+1]
        return self.xi, self.vi, self.ti

class ForwardEulerVelocity(IntegrationMethod):
    """
    Returnerer posisjon og hastighet sammen.
    """
    def __call__(self, f, x0, v0, t, steps=1):
        u = np.empty([t.size, 2] + list(x0.shape))
        u[0]      = np.array([x0, v0])
        current_u = u[0].copy()
        @jit
        def integrate(current_u):
            k = 0
            for k in range(1, t.size):
                dt = (t[k] - t[k-1]) / steps
                for i in range(steps):
                    current_u[0] = current_u[0] + dt*current_u[1]
                    current_u[1] = current_u[1] + dt*f(current_u[0], current_u[1])
                u[k] = current_u.copy()
            return u, k
        u, k = integrate(current_u)
        self.xi = u[:k+1,0,...]
        self.vi = u[:k+1,1,...]
        self.ti = t[:k+1]
        return self.xi, self.vi, self.ti

class Verlet(IntegrationMethod):
    """
    Returnerer bare posisjon.
    """
    def __call__(self, f, x0, v0, t, steps=1):
        #@jit
        def integrate(f, x0, v0, t, steps):
            x = np.empty([t.size] + list(x0.shape))
            x[0]   = x0
            prev_x = x[0]
            dt = (t[1] - t[0]) / steps
            current_x = x0 + v0*dt + 0.5*f(x[0])*dt**2
            k = 0
            for k in range(1, t.size):
                dt = (t[k] - t[k-1]) / steps
                for i in range(steps):
                    next_x = 2*current_x - prev_x + f(current_x)*dt**2
                    prev_x = current_x
                    current_x = next_x
                x[k] = current_x
            return x, k
        x, k = integrate(f, x0, v0, t, steps)
        self.xi = x[:k+1]
        self.ti = t[:k+1]
        self.vi = np.array([v0])
        return self.xi, self.vi, self.ti


class VelocityVerlet(Verlet):
    def __call__(self, f, x0, v0, t, steps=1):
        super(VelocityVerlet, self).__call__(f, x0, v0, t, steps=steps)
        self.vi = np.empty(self.xi.shape)
        for i in range(1, len(self.xi) - 1):
            self.vi[i] = (self.xi[i+1] - self.xi[i-1]) / (self.ti[i+1] - self.ti[i-1])
        self.vi[0]  = 2 * (self.xi[1]  - self.xi[0])  / (self.ti[1]  - self.ti[0])  - self.vi[1]
        self.vi[-1] = 2 * (self.xi[-1] - self.xi[-2]) / (self.ti[-1] - self.ti[-2]) - self.vi[-2]
        return self.xi, self.vi, self.ti
