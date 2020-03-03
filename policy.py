import random
import numpy as np

class MasterDockPolicy():
    def __init__(self):
        pass

    def rollout(self):
        return None

    def collect_rollout(self, buffer):
        pass

    def __call__(self, smile):
        return True

class DockPolicy():
    def __init__(self):
        pass

    def rollout(self):
        return None

    def collect_rollout(self, buffer):
        pass

    def __call__(self, smile):
        return True

class MasterMinimizePolicy:
    def __init__(self, max=100000):
        self.buffer = []

    def __call__(self, rollout):
        self.buffer.append(rollout)
        return np.mean(self.buffer)


class MinimizePolicy:
    def __init__(self, max=100000):
        self.buffer = []
        self.cutoff = None
        self.max = max
        self.bufsize = len(self.buffer)

    def rollout(self, quantile=0.1):
        return np.quantile(self.buffer, quantile)

    def collect_rollout(self, buffer):
        self.cutoff = buffer

    def __call__(self, smile, dockscore):
        if len(self.buffer) <= self.max:
            self.buffer.append(dockscore)
            self.bufsize += 1
        else:
            self.buffer[random.randint(0, self.max)] = dockscore

        if self.bufsize <= 5:
            return True
        elif self.cutoff is None:
            return dockscore <= self.rollout()
        else:
            return dockscore <= self.cutoff


class MMGBSAPolicy:
    def __init__(self):
        self.cutoff = None

    def rollout(self):
        return self.cutoff

    def collect_rollet(self, buffer):
        self.cutoff = buffer

    def __call__(self, smile, dockscore):
        if self.cutoff is not None:
            return dockscore <= self.cutoff
        else:
            return True
