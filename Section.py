# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 17:18:41 2016

@author: claudio
"""
from pint import UnitRegistry
import sympy
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


class Section:
    def __init__(self, stringers, panels):
        self.g = nx.Graph()
        self.set_stringers(stringers)
        self.set_panels(panels)
        self.cg = self.compute_cg()
        
        self.Ixx0, self.Iyy0, self.Ixy0, self.α0 = self.compute_inertia(self.cg, "ip")
        
        self.compute_princip_coords(self.α0)
        
        self.Ixx, self.Iyy, self.Ixy, self.θ = self.compute_inertia([0, 0], "pos")
        
        self.test_simmetry()
        
        #display(self.Ixx, self.Iyy, self.Ixy, self.θ)
        
    def set_stringers(self,stringers):
        for k,v in stringers.items():
            self.g.add_node(k, ip=v[0],area=v[1])
    def set_panels(self,panels):
        for k,v in panels.items():
            self.g.add_edge(k[0],k[1], thickness=v)
    def compute_cg(self):
        x_cg = sum([self.g.node[n]["area"]*self.g.node[n]["ip"][0] for n in self.g.nodes()]) / \
                    sum([self.g.node[n]["area"] for n in self.g.nodes()])
        
        y_cg = sum([self.g.node[n]["area"]*self.g.node[n]["ip"][1] for n in self.g.nodes()]) / \
                    sum([self.g.node[n]["area"] for n in self.g.nodes()])
        
        return [sympy.simplify(sympy.expand(x_cg)), sympy.simplify(sympy.expand(y_cg))]
    
    def compute_inertia(self, point, ct):
        Ixx = sum([self.g.node[n]["area"]*(self.g.node[n][ct][1]-point[1])**2 for n in self.g.nodes()])
        Iyy = sum([self.g.node[n]["area"]*(self.g.node[n][ct][0]-point[0])**2 for n in self.g.nodes()])
        Ixy = sum([self.g.node[n]["area"]*(self.g.node[n][ct][0]-point[0])* \
                                          (self.g.node[n][ct][1]-point[1]) for n in self.g.nodes()])
        Ixx = sympy.simplify(sympy.expand(Ixx))
        Iyy = sympy.simplify(sympy.expand(Iyy))
        Ixy = sympy.simplify(sympy.expand(Ixy))
        
        α = 0.5*sympy.atan(2*Ixy/(Iyy-Ixx))
        
        return [Ixx, Iyy, Ixy, α]
        
    def compute_princip_coords(self, α):
        R = sympy.Matrix([[sympy.cos(α), sympy.sin(α)],[-sympy.sin(α), sympy.cos(α)]])
        
        for n in self.g.nodes():
            b = sympy.Matrix([self.g.node[n]["ip"][0]-self.cg[0], self.g.node[n]["ip"][1]-self.cg[1]])
            self.g.node[n]["pos"] = R*b
            
            # display(self.g.node[n]["pos"])
        
        #print(Ixx, Iyy, Ixy, α)
        
        #display(Ixx, Iyy, Ixy, α)
        
    def symmetric_nodes(self, lon, coord):
        # print("Seraching among: ",lon)
        if not lon:
            return True
        else:
            for ii, nn in enumerate(lon):
                if self.g.node[nn]["pos"][coord] == 0:
                    self.symmetry[coord]["nodes"].append( (nn,nn) )                    
                    lon.pop(ii)
                    return self.symmetric_nodes(lon, coord)
                else:
                    for jj in range(ii+1,len(lon)):
                        if self.g.node[nn]["pos"][coord] == -self.g.node[lon[jj]]["pos"][coord] and \
                        self.g.node[nn]["pos"][1-coord] == self.g.node[lon[jj]]["pos"][1-coord] and \
                        self.g.node[nn]["area"] == self.g.node[lon[jj]]["area"]:
                            self.symmetry[coord]["nodes"].append( (nn,lon[jj]) )
                            lon.pop(ii)
                            lon.pop(jj-1)
                            return self.symmetric_nodes(lon, coord)
                            
        return False

    def symmetric_edges(self, loe, coord):
        
        if not loe:
            return True
        else:
            for ii, ee in enumerate(loe):
                if self.g.node[ee[0]]["pos"][coord] == -self.g.node[ee[1]]["pos"][coord] and self.g.node[ee[0]]["pos"][1-coord] == self.g.node[ee[1]]["pos"][1-coord] or \
                self.g.node[ee[0]]["pos"][coord] == 0 and self.g.node[ee[1]]["pos"][coord] == 0:
                    self.symmetry[coord]["edges"].append( (ee,ee) )
                    loe.pop(ii)
                    return self.symmetric_edges(loe, coord)
                else:
                    for jj in range(ii+1,len(loe)):
                        if (self.g.node[ee[0]]["pos"][coord] == -self.g.node[e1[0]]["pos"][0] and \
                        self.g.node[ee[0]]["pos"][1-coord] == self.g.node[e1[0]]["pos"][1] and \
                        self.g.node[ee[1]]["pos"][coord] == -self.g.node[e1[1]]["pos"][0] and \
                        self.g.node[ee[1]]["pos"][1-coord] == self.g.node[e1[1]]["pos"][1] ) or ( \
                        self.g.node[ee[0]]["pos"][coord] == -self.g.node[e1[1]]["pos"][0] and \
                        self.g.node[ee[0]]["pos"][1-coord] == self.g.node[e1[1]]["pos"][1] and \
                        self.g.node[ee[1]]["pos"][coord] == -self.g.node[e1[0]]["pos"][0] and \
                        self.g.node[ee[1]]["pos"][1-coord] == self.g.node[e1[0]]["pos"][1] ) :
                            
                    
        
    def test_simmetry(self):
        
        #loe = self.g.edges()

        # dictionary containing matched pairs, if any        
        self.symmetry = [{"nodes": [], "edges": []}, {"nodes": [], "edges": []}]
        
        if self.symmetric_nodes(self.g.nodes(),0):
            print("True for X!")
        else:
            print("False for X!")

        if self.symmetric_nodes(self.g.nodes(),1):
            print("True for Y!")
            print(self.symmetry)
        else:
            print("False for Y!")
        
        if not self.symmetry[0]["nodes"]:
            symmetric_edges(self.g.edges(),0)
            
    def plot_section(self,):
        pass
    def detect_cylces(self):
        self.cycles = nx.cycle_basis(self.g)