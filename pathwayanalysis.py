import numpy as np
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout
from graphviz import Digraph
from IPython.display import Image
import networkx as nx
import requests

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
    
    G2 = nx.DiGraph()
    
    start = data[list(data.keys())[0]]
    edge = data[list(data.keys())[1]]
    end = data[list(data.keys())[2]]
    
    unique_node = np.unique(start + end)

    for node in unique_node:
        G.node(node)
        G2.add_node(node)
    G2.nodes()
    
    for i,j in enumerate(start):
        if 'up' in edge[i]:
            G.attr('edge',color='black',style='solid',arrowhead='vee')
        else:
            G.attr('edge',color='red',style='dashed',arrowhead='tee')
        G.edge(j,end[i])
        G2.add_edge(j,end[i],label=edge[i])

    G2.edges()
    
    G.view()
    G.render('all',view=True)
    Image('all.png')

    

    return G2

def string_api(method,
               identifier,
               NoOfLim):
  '''
  methods should be among these
  1. interaction_partners
  2. network
  3. functional_annotation
  4. enrichment
  '''
  string_api_url = "https://string-db.org/api"
  output_format = "tsv"
  
  # Construct URL for STRINGdb
  request_url = "/".join([string_api_url, output_format, method])
  
  outfile = open(method+'_output-RAW.txt', 'w')

  params = {
            "identifiers" : "\r".join(["s"]), # your protein list
            "species" : 9606, # species NCBI identifier  #homo sapiens
            "echo_query" : 1, # see your input identifiers in the output  #???
            "limit" : NoOfLim # this determines the size of searching expansion - It seems a little bit arbitrary to determine the optimum number for the searching expansion 
          }
  params["identifiers"]=identifier #file to write to, named based on method used above

  # Call STRING
  results = requests.post(request_url, data=params)

  # Read and parse the results (write to outfile)
  # also, read the nodes into a dictionary similar to sif_dict: key=node, value=list of connections
  string_dict = {}
  prev_line = '' #Confirming that the line is not empty -David
  for line in results.text.strip().split("\n"):
    if line != prev_line:
      l = line.split("\t")
      outfile.write(str(l) + '\n')
      #^^ Checks the entire string_api and list (node1, node2)
      node1 = l[2]
      node2 = l[3]
      if node1 in string_dict.keys(): #if node1 already exists as key, add node2 to corresponding list
        string_dict[node1].append(node2)
      else: #if node1 doesn't exist as key, create new entry and a list containing node2 as the value
        string_dict[node1] = [node2]
    prev_line = line
  outfile.close()
  string_dict.pop('preferredName_A') #remove the header from the dictionary      

  string_output_file = open('STRING_edge_list.txt','w')
  string_set = set()

  #create edge list in a set form for STRING database retrieval
  for key in string_dict:
    # iterate through set
    for i in string_dict[key]:
      string_set.add((key, i))

  #write out to files
  for i in string_set:
    string_output_file.write(str(i[0]) + ',' + str(i[1])+ '\n')

  return string_set, results

