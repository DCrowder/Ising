import sys
from unittest import TestCase
import numpy as np
import copy
import Ising


class TestIsingModel(TestCase):
    def setUp(self, n=10):
        self.ising = Ising.IsingModel(n)
        self.ising.T = 50
        self.ising.J = 1
        self.ising.H = 0
        self.n = n

    def test_flip_state(self):
        old_box = self.ising.box
        new_box = Ising.IsingModel.flip_state(old_box.copy())  # Prevents editing existing object by copying to new var

        # print old_box, "old box value"
        # print new_box, "new box value"

        flip_count = 0
        stay_count = 0
        for i in range(self.n):
            for j in range(self.n):
                if old_box[i, j] == -1*new_box[i, j]:
                    flip_count += 1
                elif old_box[i, j] == new_box[i, j]:
                    stay_count += 1

        one_flipped = flip_count == 1
        only_one_flipped = stay_count + flip_count == self.n * self.n

        self.assertTrue(one_flipped & only_one_flipped)

    def test_check_states_not_null(self):
        has_null = False
        for i in range(self.n):
            for j in range(self.n):
                if self.ising.box[i, j] is None:
                    has_null = True
                    self.assertRaises(AssertionError)
                    break

        self.assertFalse(has_null)

    def test_check_states_have_ones(self):
        # print >> sys.stderr, self.ising.box

        has_ones = True
        for i in range(self.n):
            for j in range(self.n):
                state = self.ising.box[i, j]
                if state != 1 | state != -1:
                    has_ones = False
                    self.assertRaises(AssertionError)
                    break

        self.assertTrue(has_ones)
