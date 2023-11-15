import numpy as np

class NeuralNet (object) :
    def __init__(self, num_inputs, num_outputs, num_hidden=0, recurrent=False,
                 batch_size = 1):
        self.weights = []
        self.num_inputs = num_inputs
        self.num_outputs = num_outputs
        self.num_hidden = num_hidden
        self.recurrent = recurrent
        
        # if there is a hidden layer, we feed all inputs to hidden and all
        # hidden to outputs

        # if recurrent, we feedback hidden activations to hidden 
        # layer (or output to output layer, if there is no hidden layer)
        
        # here we create data structures for weights
        # NOTE: we consider the bias to be an extra input to each layer for
        # ease of calculation
        if num_hidden == 0 :
            n = num_inputs + 1
            if recurrent :
                n += num_outputs
            m = num_outputs
            self.weights.append(np.zeros((n,m), dtype=float))

        else :
            # input to hidden matrix first
            n1 = num_inputs + 1
            if recurrent :
                n1 += num_hidden
            m1 = num_hidden
            self.weights.append(np.zeros((n1,m1), dtype=float))

            # and now hidden to output
            n2 = num_hidden + 1
            m2 = num_outputs

            self.weights.append(np.zeros((n2, m2), dtype=float))

        self.num_params = sum([w.size for w in self.weights])
        
        # create vector for neuron activations
        self.reset(batch_size)

    def reset(self, batch_size=1) :
        self.batch_size = batch_size
        # zero out all activations
        self.activations = np.zeros((batch_size, 
                                     self.num_hidden + self.num_outputs))

    def set_params(self, params) :
        if len(params) != self.num_params :
            raise Exception("Incorrect number of params! Expected " + 
                            str(self.num_params) + ", but received " + 
                            str(len(params)))

        self.weights[0][:,:] = np.asarray(params[:self.weights[0].size]
                                          ).reshape(self.weights[0].shape)
        if self.num_hidden > 0 :
            self.weights[1][:,:] = np.asarray(params[self.weights[0].size:]
                                          ).reshape(self.weights[1].shape)

    def step(self, inputs) :
        # step the network and return the current output activations
        inputs = np.asarray(inputs, dtype=float)
        one_dimensional = (len(inputs.shape) == 1)

        # first, if just given a 1D input array, convert to a matrix
        if one_dimensional :
            inputs = inputs[None,:] # adds a second dimension
        if (inputs.shape[0] != self.batch_size) :
            raise Exception("Incorrect batch size! Should be " + 
                            str(self.batch_size) + ", but is " + 
                            str(inputs.shape[0]) ) 
        if (inputs.shape[1] != self.num_inputs) : 
            raise Exception("Incorrect number of inputs! Should be " +
                            str(self.num_inputs) + ", but is " +
                            str(inputs.shape[1]) )
        
        # second, add in biases and current hidden activations (if recurrent)
        if self.recurrent :
            if self.num_hidden > 0 :
                input_values = np.hstack( (inputs, 
                                           self.activations[:,:self.num_hidden], 
                                           np.ones((self.batch_size, 1))) )
            else :
                input_values = np.hstack( (inputs, 
                                           self.activations[:,:], 
                                           np.ones((self.batch_size, 1))) )
        else :
            input_values = np.hstack( (inputs, np.ones((self.batch_size, 1))) )
        
        # now calculate new activations by taking dot product and applying sigmoid
        new_activations = 1. / (1 + np.exp(-1.0 * np.dot(input_values, 
                                                         self.weights[0])))
        if self.num_hidden == 0 :
            self.activations[:,:] = new_activations
        else :
            # use old hidden activations to compute outputs
            hiddens = np.hstack( (self.activations[:,:self.num_hidden], 
                                  np.ones((self.batch_size,1))))
            self.activations[:,self.num_hidden:] = (
                     1. / (1 + np.exp(-1.0 * np.dot(hiddens, self.weights[1]))))
            self.activations[:,:self.num_hidden] = new_activations


        # finally, convert back to 1D, if that's what we started with
        # otherwise return the matrix
        if one_dimensional :
            return self.activations[0,self.num_hidden:]
        else :
            return self.activations[:, self.num_hidden:]
