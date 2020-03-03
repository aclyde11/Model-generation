import numpy as np

class DockPolicy():
    def __init__(self):
        pass

    def __call__(self, smile):
        return True

class MinimizePolicy():
    def __init__(self):
        self.buffer = np.zeros(10)
        self.buffer_sat = True
        self.pos = 0
        

    def __call__(self, smile, dockscore):
        v = True

        if self.buffer_sat:
            v =  dockscore <= np.quantile(self.buffer, 0.1)

        self.buffer[self.pos] = dockscore
        self.pos = (self.pos + 1) % 10
        if self.pos == 0:
            self.buffer_sat = True

        return v

def mmgbsa_ns_policy(smile, dock_score, minimize_score):
    return False
