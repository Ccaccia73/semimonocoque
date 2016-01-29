# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 17:18:41 2016

@author: claudio
"""
#from pint import UnitRegistry
import sympy
import networkx as nx
import numpy as np
#import matplotlib.pyplot as plt


class Section:
    def __init__(self, stringers, panels):        
        self.g = nx.DiGraph()
        self.set_stringers(stringers)
        self.set_panels(panels)
        
        self.g_ind = nx.Graph(self.g)        
        
        self.cg = self.compute_cg()
        
        self.Ixx0, self.Iyy0, self.Ixy0, self.α0 = self.compute_inertia(self.cg, "ip")
        
        self.compute_princip_coords(self.α0)
        
        self.Ixx, self.Iyy, self.Ixy, self.θ = self.compute_inertia([0, 0], "pos")
        
        self.detect_cycles()
        
        self.test_simmetry()
        
        #display(self.Ixx, self.Iyy, self.Ixy, self.θ)
        
    def set_stringers(self,stringers):
        for k,v in stringers.items():
            self.g.add_node(k, ip=v[0],area=v[1])
    def set_panels(self,panels):
        for k,v in panels.items():
            self.g.add_edge(k[0],k[1], thickness=v, length = sympy.sqrt( (self.g.node[k[0]]["ip"][0] - self.g.node[k[1]]["ip"][0])**2 + (self.g.node[k[0]]["ip"][1] - self.g.node[k[1]]["ip"][1])**2  ))
    def compute_cg(self):
        x_cg = sum([self.g.node[n]["area"]*self.g.node[n]["ip"][0] for n in self.g.nodes()]) / \
                    sum([self.g.node[n]["area"] for n in self.g.nodes()])
        
        y_cg = sum([self.g.node[n]["area"]*self.g.node[n]["ip"][1] for n in self.g.nodes()]) / \
                    sum([self.g.node[n]["area"] for n in self.g.nodes()])
        
        return sympy.Matrix([sympy.simplify(sympy.expand(x_cg)), sympy.simplify(sympy.expand(y_cg))])
    
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
                            # you must pop the LATTER FIRST!                            
                            lon.pop(jj)
                            lon.pop(ii)
                            
                            return self.symmetric_nodes(lon, coord)
                            
        return False

    def symmetric_edges(self, loe, coord):
        #print("Seraching among: ",loe)
        if not loe:
            return True
        else:
            for ii, ee in enumerate(loe):
                if self.g.node[ee[0]]["pos"][coord] == -self.g.node[ee[1]]["pos"][coord] and self.g.node[ee[0]]["pos"][1-coord] == self.g.node[ee[1]]["pos"][1-coord] or \
                self.g.node[ee[0]]["pos"][coord] == 0 and self.g.node[ee[1]]["pos"][coord] == 0:
                    self.symmetry[coord]["edges"].append( (ee) )
                    loe.pop(ii)
                    return self.symmetric_edges(loe, coord)
                else:
                    for jj in range(ii+1,len(loe)):
                        if self.g.edge[ee[0]][ee[1]]["thickness"] == self.g.edge[loe[jj][0]][loe[jj][1]]["thickness"] :
                            if ( self.g.node[ee[0]]["pos"][coord] == -self.g.node[loe[jj][0]]["pos"][coord] and \
                                 self.g.node[ee[0]]["pos"][1-coord] == self.g.node[loe[jj][0]]["pos"][1-coord] and \
                                 self.g.node[ee[1]]["pos"][coord] == -self.g.node[loe[jj][1]]["pos"][coord] and \
                                 self.g.node[ee[1]]["pos"][1-coord] == self.g.node[loe[jj][1]]["pos"][1-coord] ) or ( \
                                 self.g.node[ee[0]]["pos"][coord] == -self.g.node[loe[jj][1]]["pos"][coord] and \
                                 self.g.node[ee[0]]["pos"][1-coord] == self.g.node[loe[jj][1]]["pos"][1-coord] and \
                                 self.g.node[ee[1]]["pos"][coord] == -self.g.node[loe[jj][0]]["pos"][coord] and \
                                 self.g.node[ee[1]]["pos"][1-coord] == self.g.node[loe[jj][0]]["pos"][1-coord] ) :
                                     self.symmetry[coord]["edges"].append( (ee,loe[jj]) )
                                     loe.pop(jj)
                                     loe.pop(ii)
                                     
                                     return self.symmetric_edges(loe, coord)
        return False
                            
                    
        
    def test_simmetry(self):
        
        #loe = self.g.edges()
        self.ct = sympy.zeros(2,1)
        # dictionary containing matched pairs, if any        
        self.symmetry = [{"nodes": [], "edges": []}, {"nodes": [], "edges": []}]
        
        for kk in range(2):
            if self.symmetric_nodes(self.g.nodes(),kk):
                # print("Nodes True for {} !".format("X" if kk == 0 else "Y"))
                if self.symmetric_edges(self.g.edges(),kk):
                    self.ct[kk] = 0
                else:
                    self.ct[kk] = self.compute_shear_center(1-kk)
            else:
                self.ct[kk] = self.compute_shear_center(1-kk)
                #print("Edeges True for {} !".format("X" if kk == 0 else "Y"))
                #else:
                #    print("Edeges False for {} !".format("X" if kk == 0 else "Y"))
        #else:
        #    print("Nodes False for {} !".format("X" if kk == 0 else "Y"))
        
        #print(self.symmetry)
    def compute_shear_center(self, coord):
        nq = len(self.g.edges())+1
        
        self.T = sympy.zeros(nq,1)
        
        self.A = sympy.zeros(nq)

        edgedict = dict(zip(self.g.edges(), range(nq)))
        #edgedictback = dict(zip(range(nq), self.g.edges()))
        
        for rowi, nn in enumerate(self.g.nodes()[:-1]):
            for pn in self.g.predecessors(nn):
                self.A[rowi,edgedict[(pn,nn)] ] = -1
            for sn in self.g.successors(nn):
                self.A[rowi,edgedict[(nn,sn)] ] = 1
                
            self.T[rowi] = -coord/self.Ixx*self.g.node[nn]["area"]*self.g.node[nn]["pos"][1]-(1-coord)/self.Iyy*self.g.node[nn]["area"]*self.g.node[nn]["pos"][0]
         
        self.T[-(1+len(self.cycles))] = 0
        for ee in self.g.edges():
            self.A[-(1+len(self.cycles)),edgedict[ee]] = self.compute_2Omega_i(*ee, True)
        self.A[-(1+len(self.cycles)),-1] = -1

        for c_count in range(len(self.cycles)):
            cycle_nodes = self.cycles[c_count]
            cycle_nodes.append(self.cycles[c_count][0])
            for ci in range(len(cycle_nodes)-1):
                edge = (cycle_nodes[ci],cycle_nodes[ci+1])
                rev_edge = (cycle_nodes[ci+1],cycle_nodes[ci])
                if edge in self.g.edges():
                    self.A[-1-c_count,edgedict[edge]] = self.g.edge[edge[0]][edge[1]]['length'] / self.g.edge[edge[0]][edge[1]]['thickness']
                elif (rev_edge) in self.g.edges():
                    self.A[-1-c_count,edgedict[rev_edge]] = -self.g.edge[edge[1]][edge[0]]['length'] / self.g.edge[edge[1]][edge[0]]['thickness']
                else:
                    print("Problem?")

                
            
            
        
        tempq = self.A.LUsolve(self.T)
        
        return sympy.simplify(tempq[-1])
        
    def compute_Jt(self):
        pass
    
    def set_loads(self, _Tx, _Ty, _Nz, _Mx, _My, _Mz):
        self.Tx = _Tx
        self.Ty = _Ty
        self.Nz = _Nz
        self.Mx = _Mx
        self.My = _My
        self.Mz = _Mz
    
    
    def compute_stringer_actions(self):
        self.N = {}
        Atot = sum([self.g.node[n]["area"] for n in self.g.nodes()])
        for n in self.g.nodes():
            Ai = self.g.node[n]["area"]
            self.N[n] = self.Nz*Ai/Atot +self.Mx/self.Ixx*Ai*self.g.node[n]["pos"][1]-self.My/self.Iyy*Ai*self.g.node[n]["pos"][0]
    
    def compute_panel_fluxes(self):
        
        #self.q = {}
        
        nq = len(self.g.edges())
        
        self.T = sympy.zeros(nq,1)
        
        self.A = sympy.zeros(nq)
        
        edgedict = dict(zip(self.g.edges(), range(nq)))
        edgedictback = dict(zip(range(nq), self.g.edges()))
        
        for rowi, nn in enumerate(self.g.nodes()[:-1]):
            for pn in self.g.predecessors(nn):
                self.A[rowi,edgedict[(pn,nn)] ] = -1
            for sn in self.g.successors(nn):
                self.A[rowi,edgedict[(nn,sn)] ] = 1
                
            self.T[rowi] = -self.Ty/self.Ixx*self.g.node[nn]["area"]*self.g.node[nn]["pos"][1]-self.Tx/self.Iyy*self.g.node[nn]["area"]*self.g.node[nn]["pos"][0]
         
        if len(self.cycles):
            self.T[-len(self.cycles)] = self.Mz
            for ee in self.g.edges():
                self.A[-len(self.cycles),edgedict[ee]] = self.compute_2Omega_i(*ee, False)
        
        if len(self.cycles) > 1:
            c = []
            for c_count in range(len(self.cycles)-1):
                c[0] = self.cycles[c_count]
                c[0].append(self.cycles[c_count][0])
                c[1] = self.cycles[c_count+1]
                c[1].append(self.cycles[c_count+1][0])
                a_2Omega_k = [0,0]
                
                #compute area of 2 adjacent cells
                for i_om in range(2):
                    for i_node in range(len(c[i_om])-1):
                        a_2Omega_k[i_om] += self.compute_2Omega_i(c[i_om][i_node],c[i_om][i_node+1],False)

                    edge = (c[i_om][i_node],c[i_om][i_node+1])
                    rev_edge = (c[i_om][i_node+1],c[i_om][i_node])
                    if edge in self.g.edges():
                        self.A[-1-c_count,edgedict[edge]] = (-1)**i_om*self.g.edge[edge[0]][edge[1]]['length'] / self.g.edge[edge[0]][edge[1]]['thickness']
                    elif (rev_edge) in self.g.edges():
                        self.A[-1-c_count,edgedict[rev_edge]] = -(-1)**i_om*self.g.edge[edge[1]][edge[0]]['length'] / self.g.edge[edge[1]][edge[0]]['thickness']
                    else:
                        print("Problem?")


        
        tempq = self.A.LUsolve(self.T)
        
        self.q = {edgedictback[i]: sympy.simplify(tempq[i]) for i in edgedictback}
        
            
    def compute_2Omega_i(self,n1,n2, use_cg):
        v1 = sympy.zeros(3,1)
        v2 = sympy.zeros(3,1)
        
        for i in range(2):
            if use_cg:
                v1[i] = self.g.node[n1]["pos"][i]
                v2[i] = self.g.node[n2]["pos"][i]
            else:
                v1[i] = self.g.node[n1]["pos"][i] - self.ct[i]
                v2[i] = self.g.node[n2]["pos"][i] - self.ct[i]
                
                
        return v1.cross(v2)[2]
            
    
    def detect_cycles(self):
        self.cycles = nx.cycle_basis(self.g_ind)