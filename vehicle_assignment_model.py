# import solver
import gurobipy as gp #import gurobipy library in Python as gp
from gurobipy import GRB
import pandas as pd #import pandas library as pd
import numpy as np #import numpy library
import os 
print(gp.gurobi.version())

#Initialize the Gurobi model
model = gp.Model()

import ast

hl = {2,3,4,5,6,7.5,8,10,12,15,20,30,60}

#read the number of the network lines from a .txt file as an integer number (dtype='int')
L = np.loadtxt('Data/lines.txt',dtype='int') 
L = np.arange(1,L+1)
L = tuple(L)

#read the vehicle capacities for different lines from a .txt file as an integer number (dtype='int')
c = np.loadtxt('Data/capacities.txt',dtype='int') 
c = tuple(c)
c = {i:c[i-1] for i in L}
print(c)

#read the value of passenger's waiting time from a .txt file
V = np.loadtxt('Data/V.txt',dtype='float') 
print(V)

#read the cost of deploying an extra public transport vehicle from a .txt file
W = np.loadtxt('Data/W.txt',dtype='float') 
print(W)

#read the in-vehicle passenger load for the vehicles of each line to conform with the pandemic-imposed capacity limit
k = np.loadtxt('Data/lines_passengerload.txt',dtype='int') 
k = tuple(k)
k = {i:k[i-1] for i in L}
k = {i:176 for i in L}
print(k)

# Read the lost revenue per km for a passenger who is refused service
M = np.loadtxt('Data/M.txt',dtype='float') 

# Read the number of available public transport vehicles at the network level
N = np.loadtxt('Data/N.txt',dtype='int') 

#read the number of stops served by each line from a .txt file as an integer number (dtype='int')
S_number = np.loadtxt('Data/line_stops.txt',dtype='int') 
print(S_number)

#read round-trip travel times in minutes, including the layover times, for each line from a .txt file
T = np.loadtxt('Data/line_traveltimes.txt') 
T = {i:T[i-1] for i in L}
print(T)

#read the network arcs that connect the stop pairs from a .txt file that contains the arcs in a dictionary form {'arc ID': ('start_stop', 'end_stop')}. 
file = open("Data/A.txt", "r")
contents = file.read(); A = ast.literal_eval(contents); file.close()
print(A)

#read the served arcs by each line in the network from a .txt file that contains the arcs in a dictionary form. g_[l,a]=1 if line l serves arc a. Otherwise, g_[l,a]=0.
file = open("Data/g.txt", "r")
contents = file.read(); g = ast.literal_eval(contents); file.close()
print(g)

#read the maximum allowed frequency per arc (in minutes) from a .txt file
file = open("Data/f_a_max.txt", "r")
contents = file.read(); f_a_max = ast.literal_eval(contents); file.close()
print(f_a_max)

#read all origin-destination pairs in the network from a .txt file
S=np.loadtxt('Data/network_stops.txt',dtype='str')
S=tuple(S)
print('S=',S)

#read the stops of each line from a .txt file in a dictionary form {'line ID': ('stop ID', 'stop ID', '...')}
file = open("Data/line_stop_IDs.txt", "r")
contents = file.read(); Sl = ast.literal_eval(contents); file.close()
print(Sl)
test= Sl[1][0]
print('test:',test)
for s in range(1,len(Sl[1])):
    print(Sl[1][s])

#read binary parameter indicating whether an OD-pair is served by line l from a .txt file. This is in dictionary form {('line ID, 'origin stop ID', 'destination stop ID'): 1}
file = open("Data/delta.txt", "r")
contents = file.read(); delta = ast.literal_eval(contents); file.close()
print(delta)

#read binary parameter indicating whether a stop s is served by line l from a .txt file.
file = open("Data/delta_tilde.txt", "r")
contents = file.read(); delta_tilde = ast.literal_eval(contents); file.close()
print(delta_tilde)

#read origin-destination demand between OD-pairs from a .txt file
file = open("Data/Bsy.txt", "r")
contents = file.read(); Bsy = ast.literal_eval(contents); file.close()
print(Bsy)

#read the distances between the o-d pairs of each line from a .txt file.
file = open("Data/d.txt", "r")
contents = file.read(); d = ast.literal_eval(contents); file.close()
print(d)

for s in S:
    print(s)

#Initialize variable x_l denoting the number of vehicles assigned to each line l
x = model.addVars(L,vtype=gp.GRB.INTEGER, lb=1, name='x')

#initialize variable h_l denoting the time headway of line l
#h = model.addVars(L,vtype=gp.GRB.CONTINUOUS, lb=0, name='h')
h = model.addVars(L,hl,vtype=gp.GRB.BINARY, name='h')

#Initialize variable Bcapital_{l,s,y} of the expected arrival rate at station s for passenger whose destination is y and are willing to use line l as non-negative continuous variable
Bcapitallsy = model.addVars(L,S,S,vtype=gp.GRB.CONTINUOUS, lb=0, name='Bcapitallsy')

#initialize variable gamma_{l,s} denoting the vehicle load of each vehicle serving line l
gamma = model.addVars(L,S,vtype=gp.GRB.CONTINUOUS, lb=0, name='gamma')

#Initialize variable b_{l,s,y} of the hourly passenger demand between stations s and y of line l that can be served by line l while comforming to the pandemic-imposed capacity limit
blsy = model.addVars(L,S,S,vtype=gp.GRB.CONTINUOUS, lb=0, name='blsy')

#Initialize variable btilde_{l,s,y} of the hourly passenger demand between stations s and y of line l that cannot be accommodated by the public transport network due to the social distancing requirements
btildelsy = model.addVars(L,S,S,vtype=gp.GRB.CONTINUOUS, lb=0, name='btildelsy')

import math

model.addConstrs(sum(h[l,i] for i in hl) == 1 for l in L)
model.addConstrs(sum(h[l,i]*i for i in hl)*x[l] >= T[l] for l in L )
model.addConstrs(gamma[l,Sl[l][0]] == sum(blsy[l,Sl[l][0],Sl[l][y]]*(1/60)*sum(h[l,i]*i for i in hl) for y in range(1,len(Sl[l]))) for l in L)
model.addConstrs(gamma[l,Sl[l][s]] == gamma[l,Sl[l][s-1]] - sum(blsy[l,Sl[l][y],Sl[l][s]]*(1/60)*sum(h[l,i]*i for i in hl) for y in range(s,0,-1)) + sum(blsy[l,Sl[l][s],Sl[l][y]]*(1/60)*sum(h[l,i]*i for i in hl) for y in range(s,len(Sl[l]))) for l in L for s in range(1,len(Sl[l])))

model.addConstr(sum(x[l] for l in L) <= N)
#model.addConstrs( h[l]*x[l] >= T[l] for l in L )
model.addConstrs(sum(g[l,a]*(x[l]/T[l]) for l in L) <=  f_a_max[a] for a in A) 
model.addConstrs(gamma[l,s]*delta_tilde[l,s] <= k[l] for l in L for s in S)
                 

#we devide by 1/60 because h[l] is expressed in minutes
#model.addConstrs(gamma[l,Sl[l][0]] == sum(blsy[l,Sl[l][0],Sl[l][y]]*(1/60)*h[l] for y in range(1,len(Sl[l]))) for l in L)
#model.addConstrs(gamma[l,Sl[l][s]] == gamma[l,Sl[l][s-1]] - sum(blsy[l,Sl[l][y],Sl[l][s]]*(1/60)*h[l] for y in range(s,0,-1)) + sum(blsy[l,Sl[l][s],Sl[l][y]]*(1/60)*h[l] for y in range(s,len(Sl[l]))) for l in L for s in range(1,len(Sl[l])))
    
model.addConstrs(blsy[l,s,y] == delta[l,s,y]*Bcapitallsy[l,s,y]-delta[l,s,y]*btildelsy[l,s,y] for l in L for s in S for y in S)
model.addConstrs(sum(delta[l,s,y]*Bcapitallsy[l,s,y] for l in L) == Bsy[s,y] for s in S for y in S)

#Declare objective function
obj = W*sum(x[l] for l in L) + sum(sum(sum(V*blsy[l,Sl[l][s],Sl[l][y]]*(1/60)*sum(h[l,i]*i for i in hl)+M*d[l,Sl[l][s],Sl[l][y]]*btildelsy[l,Sl[l][s],Sl[l][y]] for y in range(0,len(Sl[l])) if y>s)for s in range(0,len(Sl[l])-1)) for l in L)
#obj = W*sum(x[l] for l in L) + sum(sum(sum(V*blsy[l,Sl[l][s],Sl[l][y]]*(1/60)*h[l]+M*d[l,Sl[l][s],Sl[l][y]]*btildelsy[l,Sl[l][s],Sl[l][y]] for y in range(0,len(Sl[l])) if y>s)for s in range(0,len(Sl[l])-1)) for l in L)

#Add objective function to model and declare that we solve a minimization problem
model.setObjective(obj,GRB.MINIMIZE)

model.params.NonConvex = 2 #allow to handle quadratic equality constraints - which are always non-convex
model.optimize()
print('status',GRB.Status.OPTIMAL)
model.printQuality()
for v in model.getVars():
    print('%s %g' % (v.varName, v.x))