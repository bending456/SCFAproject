import numpy as np
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout
from graphviz import Digraph
from IPython.display import Image

def dataextract(filename):
    
    rawdata = open(filename,'r')
    data = {}
    counter = 0 
    for line in rawdata: 
        if line.strip():
            if counter == 0:
                counter = counter + 1 
                line = line.strip("\n' '")
                data_name = line.split(",")
                for i,j in enumerate(data_name):
                    data[j]=[]

            else:
                line = line.strip("\n' '")
                line = line.split(",")
                key_list = list(data.keys()) 
                counter2 = 0
                for i in np.arange(len(key_list)):
                    data[key_list[i]].append(line[i])
                    counter2 += 1

    return data


def drawerForAll(data):
    G = Digraph('G', filename='all.gv',format='png',
              node_attr={'color': 'lightblue2', 'style': 'filled','shape':'box'})
    start = data[data.keys[0]]
    edge = data[data.keys[1]]
    end = data[data.keys[2]]
    
    unique_node = np.unique(start + end)
    for node in unique_node:
        G.node(node)
    for i,j in enumerate(start):
        if 'up' in edge[i]:
            G.attr('edge',color='black',style='solid',arrowhead='vee')
        else:
            G.attr('edge',color='red',style='dashed',arrowhead='tee')
        G.edge(j,end[i])
    
    G.view()
    G.render('all',view=True)
    Image('all.png')

    return

