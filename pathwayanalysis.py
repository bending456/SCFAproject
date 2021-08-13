import numpy as np

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