import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class run():
    def __init__(self, run_path):
        self.run_path = run_path
        
    def compare_params(self, params):
        print(params)
        pass
        

class chain():
    def __init__(self, chain_path):
        self.chain_path = chain_path
        
    def get_params(self):
        pass
    
    def plot_corner(self):
        pass
        
        
class transit():
    def __init__(self, chain_path, transit_idx):
        self.chain_path = chain_path
        self.transit_idx = transit_idx
        
    def plot_orbit(self):
        pass
    
    def get_transit_params(self):
        pass
    
    def plot_transit(self):
        pass