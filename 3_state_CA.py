import random
import numpy as np
import matplotlib.pyplot as plt


def random_string(length):
    """
    Returns a random bit string of the given length. 
    
    Parameters
    ----------
    length: int
        Posivite integer that specifies the desired length of the bit string.
        
    Returns
    -------
    out: list
        The random bit string given as a list, with int elements.
    """
    if not isinstance(length, int) or length < 0:
        raise ValueError("input length must be a positive ingeter")
    return [random.randint(0, 1) for _ in range(length)]


def random_Tern_string(length):
    """
    Returns a random bit string of the given length. 
    
    Parameters
    ----------
    length: int
        Posivite integer that specifies the desired length of the bit string.
        
    Returns
    -------
    out: list
        The random bit string given as a list, with int elements.
    """
    if not isinstance(length, int) or length < 0:
        raise ValueError("input length must be a positive ingeter")
    return [random.randint(0, 2) for _ in range(length)]


class ECA(object):
    """
    Elementary cellular automata simulator.
    """

    def __init__(self, rule_number, initial_condition, num_neighborhoods, num_states):
        """
        Initializes the simulator for the given rule number and initial condition.

        Parameters
        ----------
        rule_number: int
            Integer value between 0 and 255, inclusive. Specifies the ECA lookup table
            according to the Wolfram numbering scheme.
        initial_condition: list
            Ternary string used as the initial condition for the ECA. Elements of the list
            should be ints.
        num_neighborhoods: int
            The number of possible neighborhoods, for binary, we had 8, for ternary, 9.
        num_states: int
            The number of states for each cell in our CA, 2 or 3.

        Attributes
        ----------
        lookup_table: dict
            Lookup table for the ECA given as a dictionary, with neighborhood tuple keys.
        initial: array_like
            Copy of the initial conditions used to instantiate the simulator
        spacetime: array_like
            2D array (list of lists) of the spacetime field created by the simulator.
        current_configuration: array_like
            List of the spatial configuration of the ECA at the current time
        """
        # we will see a cleaner and more efficient way to do the following when we introduce numpy
        for i in initial_condition:
            if i not in [0, 1, 2]:
                raise ValueError("initial condition must be a list of 0s and 1s")

        self.num_states = num_states
        self.rule_number = rule_number
        self.num_neighborhoods = num_neighborhoods
        self.lookup_table = self.generate_lookup_table()
        self.initial = initial_condition
        self.spacetime = [initial_condition]
        self.current_configuration = initial_condition.copy()
        self._length = len(initial_condition)

    def change_base(self):
        """
        Change a number into base 2 or base 3, with 8 or 9 digits respectively.
        """
        if self.num_states == 2:
            rule_remaining = self.rule_number
            intList = []
            for i in range(self.num_neighborhoods):
                if rule_remaining / 2 ** (self.num_neighborhoods - i - 1) >= 1:
                    intList.append(1)
                    rule_remaining += -(2 ** (self.num_neighborhoods - i - 1))
                    continue
                if rule_remaining / 2 ** (self.num_neighborhoods - i - 1) < 1:
                    intList.append(0)
                    continue
            return intList

        if self.num_states == 3:
            # We have 9 neighborhoods, so we need 9 digits of ternary.
            rule_remaining = self.rule_number
            intList = []
            for i in range(self.num_neighborhoods):
                if rule_remaining / 3 ** (self.num_neighborhoods - i - 1) >= 2:
                    intList.append(2)
                    rule_remaining += -(3 ** (self.num_neighborhoods - i - 1)) * 2
                    continue
                if (
                    rule_remaining / 3 ** (self.num_neighborhoods - i - 1) >= 1
                    and rule_remaining / 3 ** (self.num_neighborhoods - i - 1) < 2
                ):
                    intList.append(1)
                    rule_remaining += -(3 ** (self.num_neighborhoods - i - 1)) * 1
                    continue
                if rule_remaining / 3 ** (self.num_neighborhoods - i - 1) < 1:
                    intList.append(0)
                    continue

            return intList

    def generate_lookup_table(self):
        """
        Take the rule number and the number of states to make a mapping from a
        string in base 2 or base 3 into a propagator that takes
        the system to the next timestep.
        """
        digits = self.change_base()
        neighborhoods = self.create_neighborhoods()

        mapping = dict()

        for i in range(len(digits)):
            # Here I want to flip the binary or ternary string string, such that
            # the 0,0,0 or 0,0 neighborhood lines up with the final digit.
            mapping.update({neighborhoods[i]: digits[len(digits) - i - 1]})
        return mapping

    def create_neighborhoods(self):
        """
        Create the allowed neighborhoods for either 2 or 3 states.
        """
        if self.num_states == 2:
            nbhds = [
                (0, 0, 0),
                (0, 0, 1),
                (0, 1, 0),
                (0, 1, 1),
                (1, 0, 0),
                (1, 0, 1),
                (1, 1, 0),
                (1, 1, 1),
            ]
        if self.num_states == 3:
            nbhds = [
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 1),
                (1, 2),
                (2, 0),
                (2, 1),
                (2, 2),
            ]
        return nbhds

    def evolve_CA(self, num_steps):
        """
        Take the initial condition and evolve it num_steps timesteps using
        the propagator defined by the
        rule number.
        """
        initial_condition = self.initial
        current_configuration = initial_condition.copy()
        propagator = self.generate_lookup_table()

        Config_Placeholder = initial_condition.copy()
        time_series = []
        time_series.append(current_configuration)
        if self.num_states == 2:
            for _ in range(num_steps):
                for x in range(len(current_configuration)):
                    xnbhd = (
                        int(current_configuration[x - 1]),
                        int(current_configuration[x]),
                        int(current_configuration[(x + 1) % self._length]),
                    )

                    Config_Placeholder[x] = int(propagator[xnbhd])
                current_configuration = Config_Placeholder.copy()
                time_series.append(current_configuration)
        if self.num_states == 3:
            for _ in range(num_steps):
                for x in range(len(current_configuration)):
                    xnbhd = (
                        int(current_configuration[x - 1]),
                        int(current_configuration[x]),
                    )

                    Config_Placeholder[x] = int(propagator[xnbhd])
                current_configuration = Config_Placeholder.copy()
                time_series.append(current_configuration)
        return time_series

    def plot_CA(self, num_steps):
        """
        Plot the result of a CA calculation.
        """
        time_series = self.evolve_CA(num_steps)
        plt.figure(figsize=(12, 12))
        # My girlfriend wanted pink, so the colormap has pink.
        plt.imshow(time_series, cmap=plt.cm.spring, interpolation="nearest")
        plt.show()
