import numpy as np
import pandas as pd
from pysimm import lmps
from sklearn.linear_model import LinearRegression

class BilinearTgFit(object):
    def __init__(self, x, y, split=None):
        self.x = x
        self.y = y
        
        if type(split) is list and len(split) == 2:
            split_lo, split_hi = split
        else:
            split_lo = split
            split_hi = split
        
        self.split_lo = split_lo
        self.split_hi = split_hi
        
        self.lofit = LinearRegression()
        self.lofit.fit(x[x<=split_lo].reshape(-1, 1), y[x<=split_lo].reshape(-1, 1))
        self.loscore = self.lofit.score(x[x<=split_lo].reshape(-1, 1), y[x<=split_lo].reshape(-1, 1))
        
        self.hifit = LinearRegression()
        self.hifit.fit(x[x>=split_hi].reshape(-1, 1), y[x>=split_hi].reshape(-1, 1))
        self.hiscore = self.hifit.score(x[x>=split_hi].reshape(-1, 1), y[x>=split_hi].reshape(-1, 1))
        
        a = np.array([[self.lofit.coef_[0][0], -1], [self.hifit.coef_[0][0], -1]])
        b = np.array([-self.lofit.intercept_[0], -self.hifit.intercept_[0]])
        self.tg, self.tg_density = np.linalg.solve(a, b)

class Tg(object):
    def __init__(self, logfiles, temperature_range, split='auto2d', split_1d_tol=1.0, split_1d_max_attempts=50, split_2d_range=None):
        self.logfiles = logfiles
        self.temperature_range = temperature_range
        self.split = split
        self.split_1d_tol = split_1d_tol
        self.split_1d_max_attempts = split_1d_max_attempts
        self.split_2d_range = split_2d_range
        if self.split_2d_range is None:
            self.split_2d_range = self.temperature_range[10:-10]
        self.y = []
        for fname in self.logfiles:
            log = lmps.LogFile(fname)
            self.y.append(1/log.data.Density.values.reshape(len(self.temperature_range), -1).mean(axis=1))
        if self.y:
            self.yerr = np.array(self.y).std(axis=0)
            self.y = np.array(self.y).mean(axis=0)
        self.x = self.temperature_range
        if self.split == 'auto2d':
            self.find_split_2d()
        elif split == 'auto1d':
            self.find_split_1d()
        self.fit = BilinearTgFit(self.x, self.y, self.split)
        
    def find_split_2d(self):
        self.scores = pd.DataFrame(columns=['split', 'tg', 'loscore', 'hiscore'])
        for s in self.split_2d_range:
            f = BilinearTgFit(self.x, self.y, split=s)
            self.scores = self.scores.append({'split': s, 'tg': f.tg, 'loscore': f.loscore, 'hiscore': f.hiscore}, ignore_index=True)
        self.split = [self.scores.loc[self.scores['loscore'].idxmax()]['split'], self.scores.loc[self.scores['hiscore'].idxmax()]['split']]
        
    def find_split_1d(self):
        split_guess = np.median(self.x)
        for att in range(self.split_1d_max_attempts):
            f = BilinearTgFit(self.x, self.y, split_guess)
            if np.abs(f.tg-split_guess) <= self.split_1d_tol:
                print(np.abs(f.tg-split_guess))
                self.split = split_guess
                return
            else:
                split_guess = f.tg
        if np.abs(f.tg-split_guess) > self.split_1d_tol:
            print('Convergence to {} not found within {} attempts'.format(self.split_1d_tol, self.split_1d_max_attempts))
