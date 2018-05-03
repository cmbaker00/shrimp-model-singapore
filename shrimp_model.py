import numpy as np
from numpy.random import binomial

class population_model:
    def __init__(self, params, max_t = np.inf):
        self.params = params
        self.max_t = max_t

    def __iter__(self):
        self.NI = self.params['NI0']
        self.NN = self.params['NN0']
        self.p = self.params['p']
        self.T = 0

        self.rI = self.params['rI0'] + self.params['xI1']*self.params['ph'] + self.params['xI2']*self.params['temp']
        self.rN = self.params['rN0'] + self.params['xN1']*self.params['ph'] + self.params['xN2']*self.params['temp']

        self.KI = self.params['KI']
        self.KN = self.params['KN']

        self.a1 = self.params['alpha_1']
        self.a2 = self.params['alpha_2']
        return self

    def __next__(self):
        if self.T >= self.max_t:
            raise StopIteration
        self.T += 1
        self.R = binomial(self.NI, self.p)
        self.D = binomial(self.NN, self.p)

        NI_change = (self.rI - self.a1*self.NN)*(self.NI - self.R) + self.KI
        NN_change = self.rN*self.NN - self.a2*self.NN*(self.NI - self.R) + self.KN

        self.NI += NI_change
        self.NN += NN_change
        if self.NI < 0:
            self.NI = 0
        if self.NN < 0:
            self.NN = 0
        return self.NI, self.NN, self.R, self.D


def run_model(model, T):
    data = [vals for vals in model]
    output = [('Invasive abundance', 'Native abundance', 'Invasive catch', 'Native catch')] + data
    return np.asarray(output)

def default_params():

    input_params = {'NI0': 250, 'NN0': 400} #Initial conditions


    input_params['p'] = .02

    input_params['rI0'] = 0.1
    input_params['xI1'] = 0.05
    input_params['xI2'] = 0.05

    input_params['rN0'] = 0.1
    input_params['xN1'] = 0
    input_params['xN2'] = 0


    input_params['KI'] = 10
    input_params['KN'] = 5

    input_params['alpha_1'] = .001
    input_params['alpha_2'] = .001

    input_params['ph'] = .02
    input_params['temp'] = -1

    return input_params

if __name__ == "__main__":
    ph_list = [1,-1,.2,.5,-.2,-.9,.5]
    temp_list = [.4,1,2,-2,0,-2,.5]
    input_params = default_params()
    max_time = 12
    site = 0
    for ph, temp in zip(ph_list, temp_list):
        site += 1
        input_params['ph'] = ph
        input_params['temp'] = temp
        model_object = population_model(input_params, max_time)
        output = run_model(model_object, max_time)

        np.savetxt("data_site_{0}.csv".format(site), output, delimiter=",",fmt='%s')
    np.savetxt("ph_data.csv", ph_list, delimiter=",",fmt='%s')
    np.savetxt("temp_data.csv", temp_list, delimiter=",",fmt='%s')
    param_array = [(key, input_params[key]) for key in input_params if key != 'ph' and key != 'temp']
    param_array = np.asarray(param_array)
    np.savetxt("true_params.csv", param_array, delimiter=",", fmt=['%s','%s'])

