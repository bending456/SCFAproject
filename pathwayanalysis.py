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

    

    return G2, G

def string_api(method,
               identifier,
               NoOfLim,
               receptorSearch):
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
            "required_score" : 0,
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

  Receptors_NCBI = dataextract('Receptor_NCBI.csv') # ['Symbol']
  Receptors_NCBI_list = Receptors_NCBI['Symbol']
  Receptor_List = ['CCR3','CD40', 'CD81', 'EGFR', 'FAS', 'IFNGR1',
                   'IFNGR2', 'IL10RA', 'IL1R1', 'IL4R', 'IL6R', 'LIFR', 'LTBR', 
                   'TLR4', 'TNFRSF13C', 'TNFRSF1A', 'TNFRSF1B', 'CXCR4', 'CXCR2', 'CCR1', 
                   'CCR5', 'P2Y1', 'P2Y12', 'P2X7', 'FFAR2', 'FFAR3', 'HADHA','FFAR4']

  #create edge list in a set form for STRING database retrieval
  for key in string_dict:
    # iterate through set
    for i in string_dict[key]:
      if receptorSearch:

        if i not in Receptors_NCBI_list:
          if i not in Receptor_List:
            string_set.add((key, i))
      else:
        string_set.add((key,i))

  #write out to files
  for i in string_set:
    string_output_file.write(str(i[0]) + ',' + str(i[1])+ '\n')

  return string_set, results

def drawer(pairset,output):
  G = Digraph('G', filename=output+'.gv',format='png',
              node_attr={'color': 'green', 'style': 'filled','shape':'box'})

  raw_node = []
  for pair in pairset:
    raw_node.append(pair[0])
    raw_node.append(pair[1])

  unique_node = np.unique(raw_node)

  for node in unique_node:
    #G.attr('node', shape='box',style='filled',color='black')
    G.node(node)

  for pair in pairset:
    #G.attr('edge',color='black',style='solid')
    G.edge(pair[0],pair[1])

  G.view()
  G.render(output,view=True)

  Image(output+'.png')

  return

def analysisByNetworkX(G,targetGene):
  path2 = nx.all_simple_paths(G,'TLR4','Immune_response')
  path = list(path2)

  PathNum = 1
  Summary = {'Path':[],
            'PathNum':[],
            'PathLength':[]}
  for i,j in enumerate(path):
    checker = 1 
    pair = set()
    if targetGene in j: # <---- choosing specific pathways
      for n,k in enumerate(j):
        if n < len(j)-1:
          edgeFeature = G.get_edge_data(k,j[n+1])
          if 'up' in edgeFeature['label']:
            feature = 1
          else:
            feature = -1
            pair.add((k,j[n+1]))
          checker = checker * feature

      Summary['Path'].append(j)
      Summary['PathNum'].append(PathNum)
      Summary['PathLength'].append(len(j))
      if checker < 0:
        print('Path #',str(PathNum))
        print(j,": Suppresing Immune Response")
        print("Length: ",len(j))
        print(pair,'\n')
        PathNum += 1 
      else:
        print('Path #',str(PathNum))
        print(j,": Promoting Immune Response")
        print("Length: ",len(j))
        PathNum += 1
        for i in pair:
          print(i[0],"------------------|",i[1])
        print('\n')

  return Summary

def findShortestPath(Summary):
  for i,j in enumerate(Summary['PathLength']):
    if i == 0:
      short = j 
    else:
      if j < short:
        short = j
        shortNum = i
    
    if i == len(Summary['PathLength'])-1:
      print('Path #',Summary['PathNum'][shortNum])
      print(Summary['Path'][shortNum])

  return shortNum

def search(targetName,midtargetName,endtargetName,numJump,cutoff):
  [string_set, difference] = string_api('network',targetName,numJump,True)
  
  listNode = []
  for pair in string_set:
    listNode.append(pair[0])
    listNode.append(pair[1])
  
  G = nx.DiGraph()
  for node in np.unique(listNode):
    G.add_node(node)
  G.nodes()
  for pair in string_set:
    G.add_edge(pair[0],pair[1])
  G.edges()

  H = G.to_undirected()
  path = nx.all_simple_paths(H,targetName,endtargetName,cutoff=cutoff)
  for i in path:
    j = 0
    if midtargetName in i:
      print(i)
      j += 1
      if j == 1:
        outcome = i
      else:
        outcome.append(i)

  return outcome


def search2(listNode,string_set,targetName,midtargetName,endtargetName,cutoff):
  
  G = nx.DiGraph()
  for node in listNode:
    G.add_node(node)
  G.nodes()
  for pair in string_set:
    G.add_edge(pair[0],pair[1])
  G.edges()
  H = G.to_undirected()
  path = nx.all_simple_paths(H,targetName,endtargetName,cutoff=cutoff)
  j = 0
  for i in path:
    if midtargetName in i:
      if j == 0:
        outcome = [list(i)]
        j += 1
      else:
        outcome.append(list(i))
      print(outcome)

  return outcome