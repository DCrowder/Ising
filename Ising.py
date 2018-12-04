import sys

import numpy as np


class IsingModel:

    def __init__(self, n, J, H, Tk):
        """Initializes system into a pair of N states.

            Takes the parameter n:

            N = nxn dimensional spin states
        """

        # self.box = np.zeros((n, n))

        # Takes a random integer from 0 to 1, multiplies by 2,
        # and subtracts one.  This gives +1 or -1 randomly distributed.
        # Returns a numpy n-dimensional array (ndarray) of size nxn
        self.n = n
        self.J = J
        self.H = 2*mu_B*H  # Convert from A/m to joules
        self.T = k_B*Tk  # Convert from K to joules
        self.beta = 1 / (k_B * T)  # Reciprocal temperature in J^-1
        self.box = 2*np.random.randint(2, size=(n, n)) - 1

    def __name__(self):
        self.name = "ising"
        return self.name

    def calc_mag(self):
        """Calculates the magnetization of a given configuration"""
        mag = sum(self.box)
        return mag

    def calc_energy(self, box):
        """Calculates the energy of a given configuration"""

        sum1 = 0

        for i in range(self.n):
            for j in range(self.n):
                # sum all nearest neighbors, impose boundary conditions of torus
                # current site spin
                s = box[i, j]
                # Neighboring spins to the right and up.  All spins are accounted for this way
                nb = box[(i+1) % n, j] + box[i, (j+1) % n]
                sum1 += nb*s

        # Sum over all the states
        sum2 = sum(box)

        # exchange interaction energy + magnetic energy
        energy = -1*self.J*sum1 - self.H*sum2
        return energy

    def mc_move(self):
        """One monte carlo move.
            Step 2: Choose a spin at random and flip it (flip_state()).
            Step 3: Compute dE = E_trial - E_old.  This is the change in energy due to a trial flip.
            Step 4: Check if dE <= 0.  If so, accept the flip
            Step 5-7: Next Check if dE > 0.  If so, see if r <= w
                where r = rand() and w = exp(-beta*dE).
                If r <= w, accept it.  If not, reject it.
                Returns nothing, but sets box as new state if successful. """

        # Step 2: Flip a state
        test_box = flip_state(self.box)

        # Step 3: Compute dE
        E_trial = self.calc_energy(test_box)
        E_old = self.calc_energy(self.box)
        dE = E_trial - E_old

        if dE <= 0:
            self.box = test_box
        else:
            r = np.random.rand()
            w = np.exp(-self.beta*dE)
            if r <= w:
                self.box = test_box

        return


    def mc_tick(self):
        """Step 8: Repeat mc_move until all spins of the system are tested.
            One sweep counts as one unit of Monte Carlo time"""
        N, T = 64, .4




    def thermalize(self, epsilon):
        """Step 9: Repeat mc_tick until thermalization occurs (equilibrium)."""


    def compute_properties(self):
        """Step 10: Compute the physical quantities of interest in n thermalized microstates.
            Do this periodically to reduce correlation between data points"""



    def compute_averages(self):
        """Step 11: Find averages of M, magnetization and E, energy of the entire system"""


    def compute_Cv(self):
        """Calculates the specific heat"""
        term1 = exp_val(self.E)
        term2 = pow(exp_val(self.E), 2)
        Cv = self.beta/self.T * (term1 - term2)

        return Cv

    def compute_chi(self):
        """Calculates the susceptibility"""
        term1 = exp_val(self.M)
        term2 = pow(exp_val(self.M), 2)
        Cv = self.beta * (term1 - term2)


""" METHODS OUTSIDE CLASS """


def flip_state(box):
    """Flips the spin of one microstate.  Takes param box (numpy ndarray) and flips a random state.
        Keep in mind this returns the new state but does not save it to the model.
        This happens at a later point in the Monte Carlo move.
        Returns box with one flipped state."""

    # Pick a random number within the indices
    i = np.random.randint(0, n-1)
    j = np.random.randint(0, n-1)

    box[i, j] *= -1  # flip from 1 to -1 or -1 to 1

    # print "index i=", i, " and index j=", j
    return box


def exp_val(A, n):
    """Find the expectation value of the given param
        A is an iterable ndarray of values, e.g. E or M
        n is the dimension of the box"""
    return sum(A)/n



"""

Constants

"""

k_B = 1.38064852e-23  # m^2*kg*s^-2*K^-1 or J/K
mu_B = 9.274e-24  # J/T or A*m^2

"""

Initial conditions

"""

# Can be user defined
J = 1  # Exchange interaction J > 0
T = 2.26918531421  # Temperature in Kelvin
H = 0  # Magnetic field

n = 4  # Total number of indices in the box
N = n*n  # total number of states in the system




