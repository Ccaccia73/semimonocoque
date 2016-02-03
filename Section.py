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
        
        α = sympy.Rational(1,2)*sympy.atan(2*Ixy/(Iyy-Ixx))
        
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
        self.is_symmetric = [False,False]        
        
        for kk in range(2):
            if self.symmetric_nodes(self.g.nodes(),kk):
                #print("Nodes True for {} !".format("X" if kk == 1 else "Y"))
                if self.symmetric_edges(self.g.edges(),kk):
                    #print("Edeges True for {} !".format("X" if kk == 1 else "Y"))
                    #self.ct[kk] = 0
                    self.is_symmetric[kk] = True
                else:
                    print("Edeges False for {} !".format("X" if kk == 1 else "Y"))
                    self.ct[kk] = self.compute_shear_center(1-kk)
            else:
                print("Nodes False for {} !".format("X" if kk == 1 else "Y"))
                print("compute SC ")
                self.ct[kk] = self.compute_shear_center(1-kk)
                #print("Edeges True for {} !".format("X" if kk == 0 else "Y"))
                #else:
                #    print("Edeges False for {} !".format("X" if kk == 0 else "Y"))
        #else:
        #    print("Nodes False for {} !".format("X" if kk == 0 else "Y"))
        
        #self.ct[0] = self.compute_shear_center(1)
        #self.ct[1] = self.compute_shear_center(0)        
        
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
        
        # Ty*xct da un contributo positivo a Mz, ma Tx*yct dà un contributo negativo a Mz, il segno dipende dalla coordinata        
        self.A[-(1+len(self.cycles)),-1] = (-1)**(coord)

        for c_count in range(len(self.cycles)):
            cycle_nodes = self.cycles[c_count]
            if cycle_nodes[0] != cycle_nodes[-1]:
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

                
            
            
        
        self.tempq = self.A.LUsolve(self.T)
        
        return sympy.simplify(self.tempq[-1])
        
    def compute_Jt(self):
        
        # tempoorarily save load values
        tmp = sympy.zeros(6,1)
        tmp[0] = self.Tx
        tmp[1] = self.Ty
        tmp[2] = self.Nz
        tmp[3] = self.Mx
        tmp[4] = self.My
        tmp[5] = self.Mz
        
        # set torque to 1
        self.Tx = sympy.Integer(0)
        self.Ty = sympy.Integer(0)
        self.Nz = sympy.Integer(0)
        self.Mx = sympy.Integer(0)
        self.My = sympy.Integer(0)
        self.Mz = sympy.Integer(1)
        
        
        # compute fluxes with torque = 1
        self.compute_panel_fluxes()
        
        # build matrix
        
        nq = len(self.g.edges())
        edgedict = dict(zip(self.g.edges(), range(nq)))


        self.A1 = sympy.zeros(nq,nq)
        
        q_prime = sympy.zeros(nq,1)
        
        for edge in self.g.edges():
            q_prime[edgedict[edge]] = self.q[edge]
            self.A1[edgedict[edge],edgedict[edge]] = self.g.edge[edge[0]][edge[1]]['length'] / self.g.edge[edge[0]][edge[1]]['thickness']
        
        self.Jt = (q_prime.T * self.A1 * q_prime)**(-1)
        
        
        self.Tx = tmp[0]
        self.Ty = tmp[1]
        self.Nz = tmp[2]
        self.Mx = tmp[3]
        self.My = tmp[4]
        self.Mz = tmp[5]
    
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
    
    def compute_panel_fluxes(self, Lcol_i = -1):
        
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
            
            if Lcol_i == -1:
                self.T[rowi] = -self.Ty/self.Ixx*self.g.node[nn]["area"]*self.g.node[nn]["pos"][1]-self.Tx/self.Iyy*self.g.node[nn]["area"]*self.g.node[nn]["pos"][0]
            else:
                self.T[rowi] = -self.L[rowi,Lcol_i]
            
        if len(self.cycles):
            if Lcol_i ==-1:
                self.T[-len(self.cycles)] = self.Mz
            else:
                self.T[-len(self.cycles)] = 0
                
            for ee in self.g.edges():
                self.A[-len(self.cycles),edgedict[ee]] = self.compute_2Omega_i(*ee, False)
        
        if len(self.cycles) > 1:
            c = [[],[]]
            for c_count in range(len(self.cycles)-1):
                c[0] = self.cycles[c_count]
                if c[0][0] != c[0][-1]:
                    c[0].append(self.cycles[c_count][0])
                c[1] = self.cycles[c_count+1]
                if c[1][0] != c[1][-1]:
                    c[1].append(self.cycles[c_count+1][0])
                               
                a_2Omega_k = [0,0]
                #compute area of 2 adjacent cells
                for i_om in range(2):
                    for i_node in range(len(c[i_om])-1):
                        a_2Omega_k[i_om] += self.compute_2Omega_i(c[i_om][i_node],c[i_om][i_node+1],False)
                    
                    for i_node in range(len(c[i_om])-1):    
                        edge = (c[i_om][i_node],c[i_om][i_node+1])
                        rev_edge = (c[i_om][i_node+1],c[i_om][i_node])
                        if edge in self.g.edges():
                            self.A[-1-c_count,edgedict[edge]] += (-1)**i_om*self.g.edge[edge[0]][edge[1]]['length'] / self.g.edge[edge[0]][edge[1]]['thickness'] / a_2Omega_k[i_om]
                        elif (rev_edge) in self.g.edges():
                            self.A[-1-c_count,edgedict[rev_edge]] += -(-1)**i_om*self.g.edge[edge[1]][edge[0]]['length'] / self.g.edge[edge[1]][edge[0]]['thickness']  / a_2Omega_k[i_om]
                        else:
                            print("Problem?")
                            print(edge)

        
        tempq = self.A.LUsolve(self.T)
        
        self.q = {edgedictback[i]: sympy.simplify(tempq[i]) for i in edgedictback}
        
        return tempq
        
            
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
        
    def compute_L(self):
        
        nnode = len(self.g.nodes())
        
        nodedict = dict(zip(self.g.nodes(), range(nnode)))
        
        self.L = sympy.zeros(nnode,nnode-3)

        aa = sympy.symbols('a0:%d'%nnode)
        
        self.expr = []
        
        if any(self.is_symmetric):
            if self.is_symmetric[0]:
                act_nodes = self.symmetry[0]['nodes']
            else:
                act_nodes = self.symmetry[1]['nodes']
                
            self.val = sympy.zeros(nnode,2)
            
            ai = 0
            a_sym = []
            a_asym = []
            for nc in act_nodes:
                if nc[0] != nc[1]:
                    self.val[nodedict[nc[0]],1] = aa[ai]
                    self.val[nodedict[nc[1]],1] = -aa[ai]
                    a_asym.append(aa[ai])                    
                    ai +=1
                
                self.val[nodedict[nc[0]],0] = aa[ai]
                self.val[nodedict[nc[1]],0] = aa[ai]
                a_sym.append(aa[ai])
                ai += 1
            
            self.expr.append(sum(self.val[:,0]))
            self.expr.append(sum([self.val[nodedict[ni],0]*self.g.node[ni]['pos'][1] for ni in self.g.nodes() ]))
            self.expr.append(sum([self.val[nodedict[ni],1]*self.g.node[ni]['pos'][0] for ni in self.g.nodes() ]))

            self.system1=sympy.zeros(len(a_sym)+1,2)
            
            for k in range(2):
                dd = sympy.poly(self.expr[k],a_sym)
                self.system1[0:len(a_sym),k] = dd.coeffs()
            
            self.sol1 = sympy.solve_linear_system(self.system1.T, *a_sym)            
            
            self.system2=sympy.zeros(len(a_asym)+1,1)
            dd = sympy.poly(self.expr[2],a_asym)
            self.system2[0:len(a_asym),0] = dd.coeffs()
            
            self.sol2 = sympy.solve_linear_system(self.system2.T, *a_asym)            
            
            free_var = [[ai for ai in a_sym if ai not in self.sol1],[ai for ai in a_asym if ai not in self.sol2]]
            
            self.free = [[],[]]
            
            for kk in range(2):        
                for snum, sym in enumerate(free_var[kk]):
                    self.free[kk].append([])
                    for snum1, sym1 in enumerate(free_var[kk]):
                        if snum==snum1:
                            self.free[kk][snum].append((sym1,1))
                        else:
                            self.free[kk][snum].append((sym1,0))
        
            self.symm_sol_list = []
            
            for ns, ls in enumerate(self.free[0]):
                self.symm_sol_list.append([(ai,self.sol1[ai].subs(ls)) for ai in self.sol1])
                for li in ls:
                    self.symm_sol_list[ns].append(li)
            
            
            for nsol, sol_list in enumerate(self.symm_sol_list):
                self.L[:,nsol] = self.val[:,0].subs(sol_list)

                
            self.asym_sol_list = []
            
            for ns, ls in enumerate(self.free[1]):
                self.asym_sol_list.append([(ai,self.sol2[ai].subs(ls)) for ai in self.sol2])
                for li in ls:
                    self.asym_sol_list[ns].append(li)

            for nsol, sol_list in enumerate(self.asym_sol_list):
                self.L[:,-1-nsol] = self.val[:,1].subs(sol_list)


        else:
            self.val = sympy.zeros(nnode,1)
            for kk in range(nnode):
                self.val[kk] = aa[kk]

            self.expr.append(sum(self.val))
            self.expr.append(sum([self.val[nodedict[ni]]*self.g.node[ni]['pos'][1] for ni in self.g.nodes() ]))
            self.expr.append(sum([self.val[nodedict[ni]]*self.g.node[ni]['pos'][0] for ni in self.g.nodes() ]))
            
            self.system1=sympy.zeros(len(aa)+1,3)

            for k in range(3):
                dd = sympy.poly(self.expr[k],aa)
                self.system1[0:len(aa),k] = dd.coeffs()
                
            self.sol1 = sympy.solve_linear_system(self.system1.T, *aa)
            
            free_var = [ai for ai in aa if ai not in self.sol1]
            
            self.free = []
            
            for snum, sym in enumerate(free_var):
                self.free.append([])
                for snum1, sym1 in enumerate(free_var):
                    if snum==snum1:
                        self.free[snum].append((sym1,1))
                    else:
                        self.free[snum].append((sym1,0))

            self.sol_list = []
            
            for ns, ls in enumerate(self.free):
                self.sol_list.append([(ai,self.sol1[ai].subs(ls)) for ai in self.sol1])
                for li in ls:
                    self.sol_list[ns].append(li)

            for nsol, sol_list in enumerate(self.sol_list):
                self.L[:,nsol] = self.val[:,0].subs(sol_list)

    def compute_H(self):
        
        self.H = sympy.zeros(len(self.g.edges()),self.L.cols)        
        
        for ii in range(self.L.cols):
            self.H[:,ii] = self.compute_panel_fluxes(ii)
            