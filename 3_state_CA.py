def random_string(length):
    '''
    Returns a random bit string of the given length. 
    
    Parameters
    ----------
    length: int
        Posivite integer that specifies the desired length of the bit string.
        
    Returns
    -------
    out: list
        The random bit string given as a list, with int elements.
    '''
    if not isinstance(length, int) or length < 0:
        raise ValueError("input length must be a positive ingeter")
    return [random.randint(0,1) for _ in range(length)]

def random_Tern_string(length):
    '''
    Returns a random bit string of the given length. 
    
    Parameters
    ----------
    length: int
        Posivite integer that specifies the desired length of the bit string.
        
    Returns
    -------
    out: list
        The random bit string given as a list, with int elements.
    '''
    if not isinstance(length, int) or length < 0:
        raise ValueError("input length must be a positive ingeter")
    return [random.randint(0,2) for _ in range(length)]

class ECA(object):
    '''
    Elementary cellular automata simulator.
    '''
    def __init__(self, rule_number, initial_condition, numNeighborhoods, numStates):
        '''
        Initializes the simulator for the given rule number and initial condition.

        Parameters
        ----------
        rule_number: int
            Integer value between 0 and 255, inclusive. Specifies the ECA lookup table
            according to the Wolfram numbering scheme.
        initial_condition: list
            Ternary string used as the initial condition for the ECA. Elements of the list
            should be ints.
        num_Neighborhoods: int
            The number of possible neighborhoods, for binary, we had 8, for ternary, 9.
        numStates: int
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
        '''
        # we will see a cleaner and more efficient way to do the following when we introduce numpy
        for i in initial_condition:
            if i not in [0,1,2]:
                raise ValueError("initial condition must be a list of 0s and 1s")


        self.numStates = numStates
        self.ruleNumber=rule_number
        self.numNeighborhoods = numNeighborhoods
        self.lookup_table = self.generate_lookup_table()
        self.initial = initial_condition
        self.spacetime = [initial_condition]
        self.current_configuration = initial_condition.copy()
        self._length = len(initial_condition)



    def change_base(self):
        '''
        Change a number into base 2 or base 3, with 8 or 9 digits respectively.
        '''
        if self.numStates==2:
            ruleRemaining=self.ruleNumber
            intList=[]
            for i in range(self.numNeighborhoods):
                if ruleRemaining/2**(self.numNeighborhoods-i-1)>=1:
                    intList.append(1)
                    ruleRemaining += -2**(self.numNeighborhoods-i-1)
                    continue
                if ruleRemaining/2**(self.numNeighborhoods-i-1)<1:
                    intList.append(0)
                    continue
            return intList

        if self.numStates==3:
            #We have 9 neighborhoods, so we need 9 digits of ternary.
            ruleRemaining=self.ruleNumber
            intList=[]
            for i in range(self.numNeighborhoods):
                if ruleRemaining/3**(self.numNeighborhoods-i-1)>=2:
                    intList.append(2)
                    ruleRemaining += -3**(self.numNeighborhoods-i-1)*2
                    continue
                if ruleRemaining/3**(self.numNeighborhoods-i-1)>=1 and 
                    ruleRemaining/3**(self.numNeighborhoods-i-1)<2:
                    intList.append(1)
                    ruleRemaining += -3**(self.numNeighborhoods-i-1)*1
                    continue
                if ruleRemaining/3**(self.numNeighborhoods-i-1)<1:
                    intList.append(0)
                    continue

            return intList

    def generate_lookup_table(self):
        '''
        Take the rule number and the number of states to make a mapping from a
        string in base 2 or base 3 into a propagator that takes
        the system to the next timestep.
        '''
        Digits=self.change_base()
        neighborhoods=self.create_Neighborhoods()


        mapping=dict()

        for i in range(len(Digits)):
            #Here I want to flip the binary or ternary string string, such that 
            #the 0,0,0 or 0,0 neighborhood lines up with the final digit.
            mapping.update({neighborhoods[i] : Digits[len(Digits)-i-1]})
        return mapping
    def create_Neighborhoods(self):
        '''
        Create the allowed neighborhoods for either 2 or 3 states.
        '''
        if self.numStates==2:
            nbhds=[(0,0,0), (0,0,1), (0,1,0), (0,1,1), (1,0,0), (1,0,1), 
                    (1,1,0), (1,1,1)]
        if self.numStates==3:
            nbhds=[(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1),
                    (2,2)]
        return nbhds


    def evolve_CA(self, numSteps):
        '''
        Take the initial condition and evolve it numSteps timesteps using
        the propagator defined by the
        rule number.
        '''
        initial_condition = self.initial
        current_Configuration = initial_condition.copy()
        propagator=self.generate_lookup_table()

        Config_Placeholder=initial_condition.copy()
        timeSeries=[]
        timeSeries.append(current_Configuration)
        if self.numStates==2:
            for i in range(numSteps):
                for x in range(len(current_Configuration)):
                    xnbhd=(int(current_Configuration[x-1]),
                            int(current_Configuration[x]),
                            int(current_Configuration[(x+1)%length]))

                    Config_Placeholder[x]=int(propagator[xnbhd])
                current_Configuration=Config_Placeholder.copy()
                timeSeries.append(current_Configuration)
        if self.numStates == 3:
            for i in range(numSteps):
                for x in range(len(current_Configuration)):
                    xnbhd=(int(current_Configuration[x-1]),
                            int(current_Configuration[x]))

                    Config_Placeholder[x]=int(propagator[xnbhd])
                current_Configuration=Config_Placeholder.copy()
                timeSeries.append(current_Configuration)
        return timeSeries

    def plot_CA(self,numSteps):
        '''
        Plot the result of a CA calculation.
        '''
        timeSeries=self.evolve_CA(numSteps)
        plt.figure(figsize=(12,12))
        #My girlfriend wanted pink, so the colormap has pink.
        plt.imshow(timeSeries, cmap=plt.cm.spring, interpolation='nearest')
        plt.show()


