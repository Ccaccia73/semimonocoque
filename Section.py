# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 17:18:41 2016

@author: claudio
"""
#from pint import UnitRegistry
import sympy
import networkx as nx
#import numpy as np
#import matplotlib.pyplot as plt


class Section:
    """
    Class used to define a semi-monocoque section
    parameters:
     - stringers: nodes of a graph defined as a dictionary: nodes numbers as key and the following attributes as list:
                  - a tuple containing x, y coordinates
                  - a parameter defining the Area of the stringer
     - panels: edges of the graph defined as a dictionary: a tuple containing start and end node as key and the following attributes:
                  - thickness of the panel
                  
    The class asumes that all parameters are sympy symbols or numbers
    """
    def __init__(self, stringers, panels):
        
        # initialize a directed graph (keeping track of orientation is useful in some calculations)
        self.g = nx.DiGraph()
        
        # define nodes of the graph
        self.set_stringers(stringers)
        
        #define edged of the graph
        self.set_panels(panels)
        
        # build an undirected graph based on the starting graph (it is useful for cycles calculations) 
        self.g_ind = nx.Graph(self.g)
        
        # compute center of gravity
        self.cg = self.compute_cg()
        
        # compute inertial properties in original reference frame
        self.Ixx0, self.Iyy0, self.Ixy0, self.α0 = self.compute_inertia(self.cg, "ip")
        
        # given the angle above and the center of gravity, the coordinates in principal reference frame are calculated
        # and saved ad a new parameter of each node ["pos"]
        self.compute_princip_coords(self.α0)
        
        # compute inertial properties in principal reference frame
        self.Ixx, self.Iyy, self.Ixy, self.θ = self.compute_inertia([0, 0], "pos")
        
        # look for all the cycles in the graph (undirected) this data is useful in multi connected sections
        self.detect_cycles()
        
        # we test if the section is symmetric wrt x or y, it is used to compute (or not) shear center and to build L matrices
        self.test_simmetry()
        
        
    def set_stringers(self,stringers):
        """
        Populate nodes of the graph
        """
        for k,v in stringers.items():
            self.g.add_node(k, ip=v[0],area=v[1])


    def set_panels(self,panels):
        """
        populate edges of the graph
        - a "lenght" parameter is calculated and added to each edge
        """
        for k,v in panels.items():
            self.g.add_edge(k[0],k[1], thickness=v, length = sympy.sqrt( (self.g.node[k[0]]["ip"][0] - self.g.node[k[1]]["ip"][0])**2 + (self.g.node[k[0]]["ip"][1] - self.g.node[k[1]]["ip"][1])**2  ))



    def compute_cg(self):
        """
        Compute center of gravity
        - return a sympy array with x and y coordinates
        """
        x_cg = sum([self.g.node[n]["area"]*self.g.node[n]["ip"][0] for n in self.g.nodes()]) / \
                    sum([self.g.node[n]["area"] for n in self.g.nodes()])
        
        y_cg = sum([self.g.node[n]["area"]*self.g.node[n]["ip"][1] for n in self.g.nodes()]) / \
                    sum([self.g.node[n]["area"] for n in self.g.nodes()])
        
        return sympy.Matrix([sympy.simplify(sympy.expand(x_cg)), sympy.simplify(sympy.expand(y_cg))])
    
    
    def compute_inertia(self, point, ct):
        """
        Compute moments of inertia
        Parameters:
        - point = center point (can be origin or c.o.g.)
        -ct = coordinate type: "ip" for coordinates in original reference frame or "pos" for coordinates w.rt. cog
        
        Return inertial properties: Ixx, Iyy, Ixy and angle
        """
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
        """
        Rototraslation of coordinate wrt angle (given as parameter) and center of gravity
        """
        
        R = sympy.Matrix([[sympy.cos(α), sympy.sin(α)],[-sympy.sin(α), sympy.cos(α)]])
        
        for n in self.g.nodes():
            b = sympy.Matrix([self.g.node[n]["ip"][0]-self.cg[0], self.g.node[n]["ip"][1]-self.cg[1]])
            self.g.node[n]["pos"] = R*b
            
        
    def symmetric_nodes(self, lon, coord):
        """
        Recursive function used to search among nodes:
        parameters:
         - lon: list of nodes: (function starts with all nodes and are progressively cancelled)
         - x coordinate or y coordinate is used for the test
         
        nodes are popped if:
         - are on axis or x1 = -x2 and y1 = y2 (or viceversa) and A1 = A2
         
        funcion is called until the list is empty (symmetric) or not (not symmetric)
         
        a list of 2 dictionaries is populated with pair symmetric of nodes (n1,n2)
         - if a node is on the axis it is repeated (n1,n1)
        """        
        
        # print("Seraching among: ",lon)
        if not lon:
            # list empty: is symmetric
            return True
        else:
            for ii, nn in enumerate(lon):
                # node on axis
                if self.g.node[nn]["pos"][coord] == 0:
                    self.symmetry[coord]["nodes"].append( (nn,nn) )                    
                    lon.pop(ii)
                    return self.symmetric_nodes(lon, coord)
                else:
                    for jj in range(ii+1,len(lon)):
                        if self.g.node[nn]["pos"][coord] == -self.g.node[lon[jj]]["pos"][coord] and \
                        self.g.node[nn]["pos"][1-coord] == self.g.node[lon[jj]]["pos"][1-coord] and \
                        self.g.node[nn]["area"] == self.g.node[lon[jj]]["area"]:
                            # symmetric nodes found
                            self.symmetry[coord]["nodes"].append( (nn,lon[jj]) )
                            # you must pop the LATTER FIRST!                            
                            lon.pop(jj)
                            lon.pop(ii)
                            
                            return self.symmetric_nodes(lon, coord)
                            
        return False

    def symmetric_edges(self, loe, coord):
        """
        Recursive function used to search among edges:

        Edges are represented as pair of nodes (n1,n2)        
        
        parameters:
         - loe: list of edges: (function starts with all edges and are progressively cancelled)
         - x coordinate or y coordinate is used for the test
         
        edges are popped if:
         - lie on axis
         - start and end node are symmetric wrt axis
         - there is another edge with same thickness and symmetric nodes wrt axis
         
        funcion is called until the list is empty (symmetric) or not (not symmetric)
         
        a list of 2 dictionaries is populated with pair symmetric of edges ((n1,n2),(n3,3))
         - if a node is on the axis it is repeated ((n1,n2),(n1,n2))
        """        
                
        #print("Seraching among: ",loe)
        if not loe:
            # empty list: symmetric
            return True
        else:
            for ii, ee in enumerate(loe):
                # node on axis or with symmetric nodes
                # TODO: is it possible to use the list of nodes?
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
        """
        function used to test the symmetry of the section:
        
        it computes symmetry of nodes and then, if they're symmetric, looks for symmetry of edges
        
        if both nodes and edges are symmetric, the corresponding coordinate of the shear center coincides with c.o.g. (i.e. = 0)
        otherwise the coordinate is calculated
        """
        
        #initialize shear center coordinates
        self.ct = sympy.zeros(2,1)
        
        # dictionary containing matched pairs, if any        
        self.symmetry = [{"nodes": [], "edges": []}, {"nodes": [], "edges": []}]
        
        # inizialize symmetry        
        self.is_symmetric = [False,False]        
        
        for kk in range(2):
            if self.symmetric_nodes(self.g.nodes(),kk):
                #print("Nodes True for {} !".format("X" if kk == 1 else "Y"))
                if self.symmetric_edges(self.g.edges(),kk):
                    #print("Edeges True for {} !".format("X" if kk == 1 else "Y"))
                    #self.ct[kk] = 0
                    self.is_symmetric[kk] = True
                else:
                    #print("Edeges False for {} !".format("X" if kk == 1 else "Y"))
                    # compute shear center coordinate
                    self.ct[kk] = self.compute_shear_center(1-kk)
            else:
                #print("Nodes False for {} !".format("X" if kk == 1 else "Y"))
                #print("compute SC ")
                self.ct[kk] = self.compute_shear_center(1-kk)
        
        #print(self.symmetry)
        
    def compute_shear_center(self, coord):
        """
        Function that computes shear center fir coordinate "coord"
        It uses an augmented version of the formula A*q = T where:
         - A : is the matrix that indicates which nodes enter or exit a node
         - q: is the vector of unknown fluxes in edges augmented with shear center position (w.r.t. C.o.G)
         - T: is the term -Ty/Jxx*Sx_i or -Tx/Jyy*S_yi for each node (with unit shear value)
        """
        
        # all edges + unknown coordinate        
        nq = len(self.g.edges())+1
        
        self.T = sympy.zeros(nq,1)
        
        self.A = sympy.zeros(nq)
        
        # mapping between pair of edges and positions in A
        # TODO: can be globally defined? 
        edgedict = dict(zip(self.g.edges(), range(nq)))
        #edgedictback = dict(zip(range(nq), self.g.edges()))
        
        for rowi, nn in enumerate(self.g.nodes()[:-1]):
            # populqate A
            for pn in self.g.predecessors(nn):
                self.A[rowi,edgedict[(pn,nn)] ] = -1
            for sn in self.g.successors(nn):
                self.A[rowi,edgedict[(nn,sn)] ] = 1
            
            # populate T
            self.T[rowi] = -coord/self.Ixx*self.g.node[nn]["area"]*self.g.node[nn]["pos"][1]-(1-coord)/self.Iyy*self.g.node[nn]["area"]*self.g.node[nn]["pos"][0]
        
        # if the section is closed add moment equation (with Mz = 0)
        self.T[-(1+len(self.cycles))] = 0
        for ee in self.g.edges():
            # compute areas w.r.t. C.o.G. 
            self.A[-(1+len(self.cycles)),edgedict[ee]] = self.compute_2Omega_i(*ee, True)
        
        # Ty*xct da un contributo positivo a Mz, ma Tx*yct dà un contributo negativo a Mz, il segno dipende dalla coordinata        
        self.A[-(1+len(self.cycles)),-1] = (-1)**(coord)
        
        # section is multi-connected: add equations for 0 rotation for each cell
        # must detect orientation of edge to assign correct sign to each term
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

                
        # solve for unknonw fluxes and coordinate
        self.tempq = self.A.LUsolve(self.T)
        
        # return just the coordinate
        return sympy.simplify(self.tempq[-1])
        
    def compute_Jt(self):
        """
        Function that computes Polar moment of inertia
         - sets unique action Mz to 1
         - compute fluxes 
         - uses the relation 1/Jt = sum q_i*q_i*l_i/t_i
        """
        # TODO: can be optimized?
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
        
        # TODO: can be globally defined?
        nq = len(self.g.edges())
        edgedict = dict(zip(self.g.edges(), range(nq)))


        self.A1 = sympy.zeros(nq,nq)
        
        q_prime = sympy.zeros(nq,1)
        
        for edge in self.g.edges():
            q_prime[edgedict[edge]] = self.q[edge]
            self.A1[edgedict[edge],edgedict[edge]] = self.g.edge[edge[0]][edge[1]]['length'] / self.g.edge[edge[0]][edge[1]]['thickness']
        
        self.Jt = (q_prime.T * self.A1 * q_prime)**(-1)
        
        # set loads back
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
        """
        function computes stresses on stringers deriving from Mx My Nz
        """
        self.N = {}
        # coimpute total area of stringers
        Atot = sum([self.g.node[n]["area"] for n in self.g.nodes()])
        for n in self.g.nodes():
            Ai = self.g.node[n]["area"]
            self.N[n] = self.Nz*Ai/Atot +self.Mx/self.Ixx*Ai*self.g.node[n]["pos"][1]-self.My/self.Iyy*Ai*self.g.node[n]["pos"][0]
    
    def compute_panel_fluxes(self, Lcol_i = -1):
        """
        function that computes fluxes on panels
        It uses the formula A*q = T where:
         - A : is the matrix that indicates which nodes enter or exit a node
         - q: is the vector of unknown fluxes in edges
         - T: is the term -Ty/Jxx*Sx_i or -Tx/Jyy*S_yi for each node
        
        Shears are supposed to be applied in shear center
        
        if the section is closed: a moment equation is added
        if the section is multi-connected: congruence equations are added
        
        the same function is used to compute columns of the H matrix:
        if the parameters is -1:
         - T is computed as indicated above
        otherwise
         - T is -L (column i)
         
        """
        #self.q = {}
        
        nq = len(self.g.edges())
        
        self.T = sympy.zeros(nq,1)
        
        self.A = sympy.zeros(nq)
        
        # TODO: can be globaly defined?
        edgedict = dict(zip(self.g.edges(), range(nq)))
        edgedictback = dict(zip(range(nq), self.g.edges()))
        
        for rowi, nn in enumerate(self.g.nodes()[:-1]):
            # populate A
            for pn in self.g.predecessors(nn):
                self.A[rowi,edgedict[(pn,nn)] ] = -1
            for sn in self.g.successors(nn):
                self.A[rowi,edgedict[(nn,sn)] ] = 1
            
            if Lcol_i == -1:
                # normal computation of fluxes
                self.T[rowi] = -self.Ty/self.Ixx*self.g.node[nn]["area"]*self.g.node[nn]["pos"][1]-self.Tx/self.Iyy*self.g.node[nn]["area"]*self.g.node[nn]["pos"][0]
            else:
                # compute the i-th column of H
                self.T[rowi] = -self.L[rowi,Lcol_i]
        
        # closed section: add moment equation
        if len(self.cycles):
            if Lcol_i ==-1:
                # normal computation: external Torque
                self.T[-len(self.cycles)] = self.Mz
            else:
                # computation of H : zero moment
                self.T[-len(self.cycles)] = 0
                
            for ee in self.g.edges():
                # compute 2Omega for each edge wrt shear center
                self.A[-len(self.cycles),edgedict[ee]] = self.compute_2Omega_i(*ee, False)
        
        # multiconnected section
        # add equations of the type theta_dot_i = theta_dot_i+1
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
        
        # compute back fluxes as a list
        self.q = {edgedictback[i]: sympy.simplify(tempq[i]) for i in edgedictback}
        
        # return fluxes (used for H computation)
        return tempq
        
            
    def compute_2Omega_i(self,n1,n2, use_cg):
        """
        function that computes 2*area of the triangle formed by the points:
        - node n1
        - node n2
        - C.o.G. if use_cg is True or Shear center
        
        the function uses cross product
        """
        
        # TODO: can be optimized, as we use only z coord of cross product?
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
        """
        function that returns a list of all (minimum cicles in the section)
        it uses an undirected version of the graph
        """
        
        # TODO: can we define the undirected graph on the fly here?
        self.cycles = nx.cycle_basis(self.g_ind)
        
    def compute_L(self):
        """
        Compute L matrix of corrective solution:
        it works differently for symmetric and non-symmetric sections
        """
        
        nnode = len(self.g.nodes())
        
        # TODO: can it be globally defined?
        nodedict = dict(zip(self.g.nodes(), range(nnode)))
        
        self.L = sympy.zeros(nnode,nnode-3)
        
        # define as many symbols as nodes
        aa = sympy.symbols('a0:%d'%nnode)
        
        self.expr = []
        
        # symmetric sections
        if any(self.is_symmetric):
            
            # select list of symmetric nodes (y-axis symmetry is privileged)
            if self.is_symmetric[0]:
                act_nodes = self.symmetry[0]['nodes']
            else:
                act_nodes = self.symmetry[1]['nodes']
            
            # assign a symbol to each node, for symmetric and antisymmetric loads
            self.val = sympy.zeros(nnode,2)
            
            ai = 0
            a_sym = []
            a_asym = []
            for nc in act_nodes:
                
                # node on axis (zero for antisymmetric, symbol fo symmetric)                
                
                if nc[0] != nc[1]:
                    # antisymmetric
                    self.val[nodedict[nc[0]],1] = aa[ai]
                    self.val[nodedict[nc[1]],1] = -aa[ai]
                    a_asym.append(aa[ai])                    
                    ai +=1
                
                #symmetric
                self.val[nodedict[nc[0]],0] = aa[ai]
                self.val[nodedict[nc[1]],0] = aa[ai]
                a_sym.append(aa[ai])
                ai += 1
            
            # write an equation for Rz = 0 (only symmetric)
            self.expr.append(sum(self.val[:,0]))
            # write and equation for Mx = 0 (only symmetric)
            self.expr.append(sum([self.val[nodedict[ni],0]*self.g.node[ni]['pos'][1] for ni in self.g.nodes() ]))
            # write an equation for My = 0 (only antisymmetric)
            self.expr.append(sum([self.val[nodedict[ni],1]*self.g.node[ni]['pos'][0] for ni in self.g.nodes() ]))
            
            # write a system of equation for Rz = 0 and Mx = 0
            self.system1=sympy.zeros(len(a_sym)+1,2)
            
            for k in range(2):
                dd = sympy.poly(self.expr[k],a_sym)
                self.system1[0:len(a_sym),k] = dd.coeffs()
            
            # solve the first underdetermined system
            self.sol1 = sympy.solve_linear_system(self.system1.T, *a_sym)            
            
            # write a system (actually only a single equation) for My = 0
            self.system2=sympy.zeros(len(a_asym)+1,1)
            dd = sympy.poly(self.expr[2],a_asym)
            self.system2[0:len(a_asym),0] = dd.coeffs()
            
            # solve second underdetermined system
            self.sol2 = sympy.solve_linear_system(self.system2.T, *a_asym)            
            # search for free variables
            free_var = [[ai for ai in a_sym if ai not in self.sol1],[ai for ai in a_asym if ai not in self.sol2]]
            
            # buld lists of combinations of free variables
            self.free = [[],[]]
            
            # for each set of equations (symmetric and antisymmetric)
            # write every possible combination of free variables (a single 1 and the rest 0)
            for kk in range(2):        
                for snum, sym in enumerate(free_var[kk]):
                    self.free[kk].append([])
                    for snum1, sym1 in enumerate(free_var[kk]):
                        if snum==snum1:
                            self.free[kk][snum].append((sym1,1))
                        else:
                            self.free[kk][snum].append((sym1,0))
            
            # solution list of symmetric component
            self.symm_sol_list = []
            
            # assign every possible set of free variables to the underdetermined system
            for ns, ls in enumerate(self.free[0]):
                self.symm_sol_list.append([(ai,self.sol1[ai].subs(ls)) for ai in self.sol1])
                for li in ls:
                    self.symm_sol_list[ns].append(li)
            
            # polulate L
            for nsol, sol_list in enumerate(self.symm_sol_list):
                self.L[:,nsol] = self.val[:,0].subs(sol_list)

            # solution list for antisymmetric component
            self.asym_sol_list = []
            
            # assign every possible set of free variables to the underdetermined system
            for ns, ls in enumerate(self.free[1]):
                self.asym_sol_list.append([(ai,self.sol2[ai].subs(ls)) for ai in self.sol2])
                for li in ls:
                    self.asym_sol_list[ns].append(li)
            
            #populate L
            for nsol, sol_list in enumerate(self.asym_sol_list):
                self.L[:,-1-nsol] = self.val[:,1].subs(sol_list)

        # section is not symmetric
        else:
            # assign variables progressively
            self.val = sympy.zeros(nnode,1)
            for kk in range(nnode):
                self.val[kk] = aa[kk]
            
            # write equations for Rz Mx My = 0
            self.expr.append(sum(self.val))
            self.expr.append(sum([self.val[nodedict[ni]]*self.g.node[ni]['pos'][1] for ni in self.g.nodes() ]))
            self.expr.append(sum([self.val[nodedict[ni]]*self.g.node[ni]['pos'][0] for ni in self.g.nodes() ]))
            
            # write a system of equations
            self.system1=sympy.zeros(len(aa)+1,3)

            for k in range(3):
                dd = sympy.poly(self.expr[k],aa)
                self.system1[0:len(aa),k] = dd.coeffs()
            
            # solve the underdetermined system
            self.sol1 = sympy.solve_linear_system(self.system1.T, *aa)
            
            # search for free variables
            free_var = [ai for ai in aa if ai not in self.sol1]
            
            self.free = []
            
            # find every possible combination of free variables (a single 1 and all zeros)
            for snum, sym in enumerate(free_var):
                self.free.append([])
                for snum1, sym1 in enumerate(free_var):
                    if snum==snum1:
                        self.free[snum].append((sym1,1))
                    else:
                        self.free[snum].append((sym1,0))
            
            # solution list
            self.sol_list = []
            # substitute each possible combination of free variables to others
            for ns, ls in enumerate(self.free):
                self.sol_list.append([(ai,self.sol1[ai].subs(ls)) for ai in self.sol1])
                for li in ls:
                    self.sol_list[ns].append(li)
            
            # polulate L
            for nsol, sol_list in enumerate(self.sol_list):
                self.L[:,nsol] = self.val[:,0].subs(sol_list)

    def compute_H(self):
        """
        compute matrix H of corrective solution:
        - uses panel fluxes calculation for each column
        """
        
        self.H = sympy.zeros(len(self.g.edges()),self.L.cols)        
        
        for ii in range(self.L.cols):
            self.H[:,ii] = self.compute_panel_fluxes(ii)
     
    def compute_KM(self, A0, l0, t0):
        """
        compute Ktilde and Mtilde matrices in function of:
        - reference area A0
        - reference length l0
        - reference thickness t0
        
        Ktilde = L^T * [A0/A] * L
        Mtilde = H^T * [l/l0*t0/t] * H
        
        the eigenvalues of Ktilde*beta2 - Mtilde are also computed
        """
        A0_A = sympy.eye(len(self.g.nodes()))
        
        
        nq = len(self.g.edges())
        # TODO: can be globally defined?        
        edgedict = dict(zip(self.g.edges(), range(nq)))

        lt = sympy.eye(nq)
        
        for ni,node in enumerate(self.g.nodes()):
            A0_A[ni,ni] = A0/self.g.node[node]["area"]
            
        for edge in self.g.edges():
            lt[edgedict[edge],edgedict[edge]] = self.g.edge[edge[0]][edge[1]]["length"]/self.g.edge[edge[0]][edge[1]]["thickness"]*t0/l0
            
        self.Ktilde = self.L.T*A0_A*self.L
        
        self.Mtilde = self.H.T*lt*self.H
        
        Mat1 = self.Ktilde.inv()*self.Mtilde        
        
        self.β2 = Mat1.eigenvals()
        
        
        