import numpy as np
import matplotlib.pyplot as plt
from abc import *
import time

class StatModel(ABC):
    def __init__(self, N, J, h):
        self.N = N
        self.J = J
        self.h = h

    @abstractmethod
    def hamiltonian(self, state):
        raise NotImplementedError()

    @abstractmethod
    def update(self, state):
        raise NotImplementedError()

    @abstractmethod
    def initializeState(self):
        raise NotImplementedError()

    def accept_prob(self, curr_obj, new_obj, temperature):
        tmp = np.exp(-(new_obj-curr_obj)/temperature)
        return min([tmp,1.])

    def scheduleLogarithm(self, step, beta0=1.0):
        return 1/(beta0 * np.log(2+step))


    def simulated_annealing(self, maxStep, T_init=1e3, cooling_rate=0.995, seed=24, monte_carlo_step=100000):
        np.random.seed(seed)
        t = T_init
        hams = []
        ts = []
        state = self.initializeState()
        startTime = time.time()
        times = []

        for _ in range(maxStep):
            for _ in range(monte_carlo_step):
                curr_obj = self.hamiltonian(state)
                hams.append(curr_obj)
                times.append(time.time()-startTime)
                new_state = self.update(state.copy())
                new_obj = self.hamiltonian(new_state)
                if(np.random.rand() < self.accept_prob(curr_obj,new_obj,t)):
                    state = new_state

            t *= cooling_rate    # cooling


        hams.append(curr_obj);
        times.append(time.time()-startTime)
        return state, hams, times

class IsingModel(StatModel):
    def __init__(self, N, J, h):
        super().__init__(N, J, h.reshape((N,)));

    # override
    def hamiltonian(self, state):
        tmp1 = -0.5 * state.reshape((1,self.N)) @ self.J @ state.reshape((self.N,1))
        tmp2 = -state.reshape((1,self.N)) @ self.h
        return (tmp1 + tmp2)[0,0]

    # override
    def update(self, state):
        idx = np.random.randint(low=0,high=self.N)
        state[idx] = -state[idx]    # flip
        return state

    # override
    def initializeState(self):
        return np.sign(np.random.randn(self.N,))

class PottsModel (StatModel):
    def __init__(self, N, J, h, q):
        super().__init__(N, J, h.reshape((N,)))
        self.q = q

    # override
    def hamiltonian(self, state):
        # 位相表記に変える
        phase = 2*np.pi*state/self.q
        rv = phase.reshape((1,self.N))
        cv = phase.reshape((self.N,1))
        tmp1 = self.J * np.cos(cv-rv)
        tmp2 = self.h @ phase.reshape(self.N,)
        return -0.5 * tmp1.sum() - tmp2

    # override
    def update(self, state):
        idx = np.random.randint(low=0, high=self.N)
        state[idx] = np.random.randint(low=0, high=self.q)
        return state

    #override
    def initializeState(self):
        return np.random.randint(low=0, high=self.q, size=(self.N,))

class XYModel(StatModel):
    def __init__ (self, N, J, h):
        super().__init__(N, J, h.reshape((2,N)))

    # override
    def hamiltonian(self, state):
        tmp1 = -0.5 * self.J * (state.T @ state)
        tmp2 = -self.h * state
        return tmp1.sum() + tmp2.sum()

    #override
    def update(self, state):
        idx = np.random.randint(low=0, high=self.N)
        tmp = 2*np.pi*np.random.rand()
        state[0,idx] = np.cos(tmp)
        state[1,idx] = np.sin(tmp)
        #print (np.arctan2(state[1,idx], state[0,idx]))
        return state

    # override
    def initializeState(self):
        tmp = 2*np.pi*np.random.rand(1,self.N)
        return np.r_[np.cos(tmp), np.sin(tmp)]

