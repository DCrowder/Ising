import sys

import numpy as np


class IsingModel:

    def __init__(self, n):
        """Initializes system into a pair of N states.

            Takes the parameter n:

            N = nxn dimensional spin states
        """

        # self.box = np.zeros((n, n))

        # Takes a random integer from 0 to 1, multiplies by 2,
        # and subtracts one.  This gives +1 or -1 randomly distributed.
        # Returns a numpy n-dimensional array (ndarray) of size nxn
        self.box = 2*np.random.randint(2, size=(n, n)) - 1

    def __name__(self):
        self.name = "ising"
        return self.name

    def calc_mag(self):
        """Calculates the magnetization"""

    def calc_energy(self):
        """Calculates the energy"""

    @staticmethod
    def flip_state(box):
        """Flips the spin of one microstate.  Takes param box (numpy ndarray) and flips a random state.
            Keep in mind this returns the new state but does not save it to the model.
            This happens at a later point in the Monte Carlo move."""

        # Pick a random number within the indices
        i = np.random.randint(0, n-1)
        j = np.random.randint(0, n-1)

        box[i, j] = -1*box[i, j]  # flip from 1 to -1 or -1 to 1

        # print "index i=", i, " and index j=", j
        return box


"""

Constants

"""

k_B = 1.38064852e-23  # m^2*kg*s^-2*K^-1 or J/K

"""

Initial conditions

"""

T = 2.26918531421  # Temperature in Kelvin

n = 4  # Total number of indices in the box

N = n*n  # total number of states in the system

beta = 1/(k_B * T)  # Reciprocal temperature in J^-1

