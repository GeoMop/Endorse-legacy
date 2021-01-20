import numpy as np
import attr
import matplotlib.pyplot as plt


@attr.s(auto_attribs=True)
class Transport:
    velocity : float = 10
    diffusion : float = 1e-1  # Diffusion coeeficient
    rho :float = 2700  # rock density
    kd : float = 0  # distribution coefficient of linear sorption
    porosity : float = 1  # porosity
    retardation_ : float = None

    lam = 0  # decay rate

    length :float = 1000  # transport path length
    end_time : float = 250  # simulation period
    dt : float = 0.1  # number of time step (higher is better but slower)
    
    @property
    def retardation(self):
        if self.retardation_ is None:
            self.retardation_ = (1 + self.kd * self.rho * (1 - self.porosity) / self.porosity)          
        return self.retardation_
    
    @property
    def a(self):
        return self.retardation / (2 * np.sqrt(self.diffusion * self.retardation))

    @property
    def b(self):
        return self.velocity / (2 * np.sqrt(self.diffusion * self.retardation))

    @property
    def c(self):
        return self.velocity / self.diffusion

    def n_steps(self):
        return int(np.ceil(self.end_time / self.dt))
    
    def dt_(self):
        return self.end_time / self.n_steps()
    
    def times(self):
        return np.linspace(0, self.end_time, self.n_steps() + 1)
    
    
    def kernel(self, x = None):
        if x is None:
            x = self.length
        #assert(self.end_time > self.a * self.length / self.b,
        #        'Too short simulation interval.')
        
        t = self.times()
        sqrt_t = np.sqrt(t)
        ax = self.a * x
        bt = self.b * t
        e_1 = -np.square(ax - bt) / t
        e_2 = self.c * x - np.square(ax + bt) / t
        print(e_1 - e_2)
        
        normalize = 1/(2*np.sqrt(np.pi))
        kern = normalize * (
                np.exp(e_1) * (ax / t / sqrt_t + self.b / sqrt_t)
                + np.exp(e_2) * (ax / t / sqrt_t - self.b / sqrt_t)
                )
        #kern = normalize * np.exp(e_1) * ax / t / sqrt_t
        kern[0] = 0
        return kern
        
    def conc_out(self, bc_in):
        N = self.n_steps() + 1
        return np.convolve(bc_in, self.kernel(), mode='full')[:N] * self.dt_() 
    
    
trans = Transport()

print(trans.a, trans.b, trans.c)
t = trans.times()
kernel = trans.kernel()
#bc_in = np.ones_like(t)
bc_in = np.where(t < 10, 1, 0)
conc_out = trans.conc_out(bc_in)

fig, ax = plt.subplots()
ax.plot(t, bc_in, 'r')
ax.plot(t, kernel, 'b')
ax.plot(t, conc_out , 'g')

plt.show()


