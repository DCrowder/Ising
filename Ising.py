import sys

import numpy as np


class IsingModel:

    def __init__(self, n, J, H):
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
        self.box = 2*np.random.randint(2, size=(n, n)) - 1

    def __name__(self):
        self.name = "ising"
        return self.name

    def calc_mag(self):
        """Calculates the magnetization of a given configuration"""
        mag = np.sum(self.box)
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
                nb = box[(i+1) % n, j] + box[i, (j+1) % n] + box[(i-1) % n, j] + box[i, (j-1) % n]
                sum1 += nb*s

        sum1 /= 4

        # Sum over all the states
        sum2 = np.sum(box)

        # exchange interaction energy + magnetic energy
        energy = -1*self.J*sum1 - self.H*sum2
        return energy

    def mc_move(self, beta):
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
            w = np.exp(-beta*dE)
            if r <= w:
                self.box = test_box

        return


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

    return box


"""

Constants

"""

k_B = 1.38064852e-23  # m^2*kg*s^-2*K^-1 or J/K
mu_B = 9.274e-24  # J/T or A*m^2

"""

Initial conditions

"""

# Physical parameters
J = 1  # Exchange interaction J > 0
Tc = 2.26918531421  # Critical Temperature in Kelvin
H = 0  # Magnetic field

#
nt = 88 # Number of temperature points
n = 4  # Total number of indices in the box
N = n*n  # total number of states in the system
eqSteps = 1024  # Number of MC sweeps for equilibration
mcSteps = 1024  # Number of MC sweeps for calculation

# Generates nt samples from start to stop temperatures.  Returns an ndarray
T = np.linspace(1.53, 3.28, nt)

# Properties of interest
E, M, C, X = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt) # Returns ndarrays of nt samples

# Divide by number of samples and system size to get intensive values
n1, n2 = 1.0/(mcSteps*N), 1.0/(mcSteps*mcSteps*N)

# Iterate over the temperature samples
for tt in range(nt):
    E1 = M1 = E2 = M2 = 0  # Initialize averages
    config = IsingModel(n, J, H)  # Initialize states
    iT=1.0/T[tt]; iT2=iT*iT  # beta value and beta squared

    for i in range(eqSteps):    # iterate until equilibrium
        config.mc_move(iT)

    for i in range(mcSteps):
        config.mc_move(iT)

        Ene = config.calc_energy(config.box)  # calculate the energy
        Mag = config.calc_mag()   # calculate magnetization

        E1 += Ene
        M1 += Mag
        M2 += Mag*Mag
        E2 += Ene*Ene

    E[tt] = n1*E1   # Energy ndarray
    M[tt] = n1*M1   # Magnetization ndarray
    C[tt] = iT2*(n1*E2 - n2*E1*E1)
    X[tt] = iT*(n1*M2 - n2*M1*M1)

    print >> sys.stdout, 1/iT, E[tt], M[tt], C[tt], X[tt]




