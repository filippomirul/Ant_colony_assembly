from inspyred import ec
from inspyred.benchmarks import Benchmark
from utils.utils_08.ann import NeuralNet
import numpy as np

MAX_WEIGHT = 8

INPUTS = np.zeros((4,2), dtype=float)
for i in range(4) :
    INPUTS[i,:] = np.array([1 if i<2 else 0, 1 if i%2==0 else 0])

class NeuralNetworkBenchmark(Benchmark):
    """Defines the base class for Neural Network Benchmark Problems.  Other 
        neural net benchmarks should inherit from this
    """
    
    def __init__(self, net):
        self.net = net
        super(NeuralNetworkBenchmark, self).__init__(self.net.num_params)
        self.bounder = ec.Bounder([-MAX_WEIGHT] * self.dimensions, 
                                  [MAX_WEIGHT] * self.dimensions)
    
    def evaluator(self, candidates,args) :  
        return np.array(list(map(self.evaluate_single, candidates)))

    def generator(self, random, args):
        return np.asarray([random.uniform(-MAX_WEIGHT, MAX_WEIGHT) 
                    for _ in range(self.net.num_params)])       

class BaseLogicBenchmark(NeuralNetworkBenchmark):
    """Defines the base class for single time step logic 
        Neural Network Benchmark Problems. Other logic neural net benchmarks 
        should inherit from this
    """
    
    def __init__(self, net, logic_fn) :
        super(BaseLogicBenchmark, self).__init__(net)
        self.maximize = False
        self.targets = np.asarray(logic_fn(INPUTS[:,0],INPUTS[:,1]), 
                                  dtype=float)[:,None]
        
    def evaluate_single(self, candidate) :
        self.net.set_params(candidate)
        self.net.reset(4)
        # step once if no hidden
        outputs = self.net.step(INPUTS)
        # if hidden, step a second time so inputs reach outputs
        if self.net.num_hidden > 0 :
            outputs = self.net.step(INPUTS)

        # feed the inputs to the net and calculate error
        return np.sum((outputs - self.targets) ** 2)

class Or(BaseLogicBenchmark):
    """Defines Or Benchmark function, using a neural net
    """

    def __init__(self, num_hidden=0, recurrent=False) :
        super(Or, self).__init__(NeuralNet(2,1,num_hidden,recurrent),
                                  np.logical_or)
 
class And(BaseLogicBenchmark):
    """Defines OR Benchmark function, using a neural net
    """

    def __init__(self, num_hidden=0, recurrent=False) :
        super(And, self).__init__(NeuralNet(2,1,num_hidden,recurrent),
                                  np.logical_and)

class Xor(BaseLogicBenchmark):
    """Defines Xor Benchmark function, using a neural net
    """

    def __init__(self, num_hidden=0, recurrent=False) :
        super(Xor, self).__init__(NeuralNet(2,1,num_hidden,recurrent),
                                  np.logical_xor)

class TemporalLogicBenchmark(NeuralNetworkBenchmark):
    """Defines the base class for temporal logic 
        Neural Network Benchmark Problems. Other temporal logic neural net 
        benchmarks should inherit from this
    """
        
    def __init__(self, net, logic_fn) :
        super(TemporalLogicBenchmark, self).__init__(net)
        self.maximize = False
        self.targets = np.asarray(logic_fn(INPUTS[:,0],INPUTS[:,1]), 
                                  dtype=float)[:,None]
        
    def evaluate_single(self, candidate) :
        self.net.set_params(candidate)
        self.net.reset(4)
        
        # for temporal xor we care about the output after seeing both inputs
        for i in range(2) :
            # step once if no hidden
            outputs = self.net.step(INPUTS[:,i:i+1])
            # if hidden, step a second time for each input so inputs 
            # reach outputs
            if self.net.num_hidden > 0 :
                outputs = self.net.step(INPUTS[:,i:i+1])

        # feed the inputs to the net and calculate error
        return np.sum((outputs - self.targets) ** 2)  
        
class TemporalOr(TemporalLogicBenchmark):
    """Defines Temporal Or Benchmark function, using a neural net
    """

    def __init__(self, num_hidden=0, recurrent=False) :
        super(TemporalOr, self).__init__(NeuralNet(1,1,num_hidden,recurrent),
                                          np.logical_or)

class TemporalAnd(TemporalLogicBenchmark):
    """Defines OR Benchmark function, using a neural net
    """

    def __init__(self, num_hidden=0, recurrent=False) :
        super(TemporalAnd, self).__init__(NeuralNet(1,1,num_hidden,recurrent),
                                  np.logical_and)           
            

class TemporalXor(TemporalLogicBenchmark):
    """Defines Temporal Xor Benchmark function, using a neural net
    """

    def __init__(self, num_hidden=0, recurrent=False) :
        super(TemporalXor, self).__init__(NeuralNet(1,1,num_hidden,recurrent),
                                          np.logical_xor)
