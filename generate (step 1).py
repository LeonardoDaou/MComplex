##############Protein Complex identification algorithm based on DPPN########################
import sys
from math import sqrt
import pandas as pd
from optparse import OptionParser
from scipy import *
import scipy as sp

from scipy.linalg import *

from numpy import *
import numpy as np

from numpy.linalg import *
import string
import os

protein_out_file="protein.temp"
cliques_file="cliques"
ppi_pair_file="ppi.pair"
ppi_matrix_file="ppi.matrix"
new_PPI_file="new_ppi.txt"

def f_key(a):
    return(a[-1])
def getOptions():
    optparser = OptionParser(usage="%prog [options]\n-h for help")
    optparser.add_option("-p", "--PPI_file", dest="PPI_file", help="Input file of PPI data",default="krogan.txt")
    
    optparser.add_option("-e", "--Gene_expression_file", dest="EXPRESSION_file", help="Input file of gene expression data",default="series_matrix.txt")
    
    optparser.add_option("-o", "--output", dest="output", help="Output file for writing the identified complexes",default="result.txt",)
   
    optparser.add_option("-c", "--core_thes_parameter", dest="Core_thes_parameter", help="expand_thes_parameter: range from 0 to 1", default="0.1")
  
    (options, args) = optparser.parse_args()
    if not options.PPI_file:
        optparser.error("No input PPI data file")
    if not options.EXPRESSION_file:
        optparser.error("No input Gene expression data file")
    if not options.output:
        optparser.error("No output file defined")
    return options, args

if __name__ == "__main__":
    
    options, args = getOptions()

    Dic_map={}
    index=0
    All_node=set([])
    All_node_index=set([])
    go_set_list=[]
    expression_list=[]
    expression_PPI_list=[]
    All_node_go=set([])

    protein_time_list=[]
    Seed_edge_list=[]
    Cliques_list=[]
    Core_list=[]
    Complex_list=[]
    Core_protein_set=set([])
    
    Core_edge_set=set([])
    Protein_SD_weight={}
   
    Node_PPI_1=[]
    Node_PPI_2=[]

    neighbor_PPI_list=[]
    Time_protein_dic=[]

    Time_protein_window_dic=[]
    Time_num=12
    Times_for_SD=3
    
    SD3_weight=0.99
    SD2_weight=0.95
    SD1_weight=0.68
    Time_window=2 
    Half_time_window=Time_window/2
    
 
    Core_weight_thresh=float(options.Core_thes_parameter)
  
    ###########input high-throughput PPI data############
    f = open(options.PPI_file,"r")
    f_protein_out=open(protein_out_file,"w")
    
    for line in f:
        line = line.strip().split()
        if len(line)==2:
           if line[0] not in Dic_map:
                
                Dic_map[line[0]]=index
                f_protein_out.write(line[0]+"\n")
                
                
                go_set_list.append(set([]))
                expression_list.append([])
                expression_PPI_list.append([])
                
                neighbor_PPI_list.append([])
                for i in range(0,Time_num//Time_window):
                    neighbor_PPI_list[index].append(set([]))
                    
                    
                protein_time_list.append(set([])) 
                index+=1
                
           if line[1] not in Dic_map:                
                Dic_map[line[1]]=index
                f_protein_out.write(line[1]+"\n")
                
             
                go_set_list.append(set([]))
                expression_list.append([])
                expression_PPI_list.append([])
                
                neighbor_PPI_list.append([])
                for i in range(0,Time_num//Time_window):
                    neighbor_PPI_list[index].append(set([]))
                    
                protein_time_list.append(set([]))

                index+=1
                
    Protein_num=index
    f.close()
    f_protein_out.close()

    ###########input Gene expression data##################
   
    f=open(options.EXPRESSION_file,"r")
    
    
    for line in f:
        line = line.strip().split()
        if len(line)==38:
            if line [1] in Dic_map:
                
                for i in range(2,Time_num+2):
                                        
                    expression_list[Dic_map[line[1]]].append((float(line[i])+float(line[i+12])+float(line[i+24]))/3)

                for j in range(0,3):

                    for i in range(2,Time_num+2):
                                        
                        expression_PPI_list[Dic_map[line[1]]].append((float(line[i])+float(line[i+12])+float(line[i+24]))/3)
         
    f.close()


    ###############compute active time attribute for proteins###################
    for i in range(0,Time_num+1):
        Time_protein_dic.append({})

    
    for instance in Dic_map:
        
        if len(expression_list[Dic_map[instance]])>=12:
            Temp_mean_value= 0.
            Temp_sd_value=0.
            Temp_thresh_3SD=0.
            Temp_thresh_2SD=0.
            Temp_thresh_1SD=0.
            
            for j in range(0,Time_num):
                Temp_mean_value+=expression_list[Dic_map[instance]][j]
            Temp_mean_value/=Time_num

            for j in range(0,Time_num):
                Temp_sd_value+=(expression_list[Dic_map[instance]][j]-Temp_mean_value)**2
            Temp_sd_value/=(Time_num-1)
            
            Temp_thresh_3SD=Temp_mean_value+3*(Temp_sd_value**0.5)*(Temp_sd_value/(1+Temp_sd_value))
            Temp_thresh_2SD=Temp_mean_value+2*(Temp_sd_value**0.5)*(Temp_sd_value/(1+Temp_sd_value))
            Temp_thresh_1SD=Temp_mean_value+1*(Temp_sd_value**0.5)*(Temp_sd_value/(1+Temp_sd_value))
            
            for j in range(0,Time_num):
                if expression_list[Dic_map[instance]][j]>=Temp_thresh_3SD:
                    Time_protein_dic[j][instance]=SD3_weight
                   
                elif Temp_thresh_3SD>expression_list[Dic_map[instance]][j]>=Temp_thresh_2SD:
                    Time_protein_dic[j][instance]=SD2_weight
                    
                elif Temp_thresh_2SD>expression_list[Dic_map[instance]][j]>=Temp_thresh_1SD:
                    Time_protein_dic[j][instance]=SD1_weight
    for instance in Dic_map:
        if instance in Time_protein_dic[0]:
            Time_protein_dic[Time_num][instance]=Time_protein_dic[0][instance]
            
              
    
    ###########input PPI data############

    for i in range(0,Time_num//Time_window):
        Node_PPI_1.append([])
        Node_PPI_2.append([])
        Time_protein_window_dic.append({})
  
    f = open(options.PPI_file,"r")
    f_out=open(new_PPI_file,"w")
    
    for line in f:
        line = line.strip().split()
        if len(line)==2:
            if line[0] !=line[1]:
                for i in range(0,Time_num//Time_window):
                    temp_weight_0=0.
                    temp_weight_1=0.

                    for j in range(i*Time_window,i*Time_window+Time_window+1):
                        #print(j)
                        if line[0] in Time_protein_dic[j]:
                            temp_weight_0+=Time_protein_dic[j][line[0]]
                        if line[1] in Time_protein_dic[j]:
                            temp_weight_1+=Time_protein_dic[j][line[1]]
                    if temp_weight_0 !=0. and temp_weight_1 !=0.:
                        
                        Node_PPI_1[i].append(line[0])
                        
                        Node_PPI_2[i].append(line[1])
                        Time_protein_window_dic[i][line[0]]=temp_weight_0/(Time_window+1)
                        Time_protein_window_dic[i][line[1]]=temp_weight_1/(Time_window+1)
                        neighbor_PPI_list[Dic_map[line[0]]][i].add(line[1])
                        neighbor_PPI_list[Dic_map[line[1]]][i].add(line[0])

    f_out.close()
                
  
    f.close()

    ######dic_map to map_dic &&&& all_protein_set###########
    df=pd.DataFrame(Time_protein_window_dic)
    df.to_excel("win_dic.xlsx")

    Map_dic={}
    for key in Dic_map.keys():
        Map_dic[Dic_map[key]]=key
        All_node_index.add(Dic_map[key])

    with open("myfile.txt", 'w') as f: 
        for key, value in Dic_map.items(): 
            f.write('%s:%s\n' % (key, value))

    ###########a cycle divided into some period##################
  
    for Time_T in range(0,Time_num//Time_window):


        ######bulid Adj_matrix at Time_T ###########
        
        Adj_Matrix = mat(zeros((Protein_num,Protein_num), dtype = int))
        
        Adj_Matrix_Weight = mat(zeros((Protein_num,Protein_num), dtype = float))
        

        if len(Node_PPI_1[Time_T])==len(Node_PPI_2[Time_T]):
            
            for i in range(len(Node_PPI_1[Time_T])):
                if Node_PPI_1[Time_T][i] in Dic_map and Node_PPI_1[Time_T][i] in Dic_map:
                    Adj_Matrix[Dic_map[Node_PPI_1[Time_T][i]],Dic_map[Node_PPI_2[Time_T][i]]]=1
                    Adj_Matrix[Dic_map[Node_PPI_2[Time_T][i]],Dic_map[Node_PPI_1[Time_T][i]]]=1

        ################compute weight for PPIN at Time_T ################

        Weight_PPI =0.
        Weight_GO =0.
        Weight_Expression =0.
        Temp_Weight_Edge=0.
        
        threshold_r=0.5
       
        print("compute weight for PPIN")
        temp_num=0
                               
        for i in range(0,len(Node_PPI_1[Time_T])):

            if len(expression_PPI_list[Dic_map[Node_PPI_1[Time_T][i]]])>=36 and len(expression_PPI_list[Dic_map[Node_PPI_2[Time_T][i]]])>=36:
                temp_PPI_expression_avg_1=0.
                temp_PPI_expression_avg_2=0.
                temp_sum_1=0.
                temp_sum_2=0.
                temp_sum_3=0.
                temp_r=0.
                
                for j in range(Time_T*Time_window,Time_T*Time_window+Time_window+1):
                    temp_PPI_expression_avg_1+=expression_PPI_list[Dic_map[Node_PPI_1[Time_T][i]]][j]
                    temp_PPI_expression_avg_2+=expression_PPI_list[Dic_map[Node_PPI_2[Time_T][i]]][j]
                   
                temp_PPI_expression_avg_1/=Time_window
                temp_PPI_expression_avg_2/=Time_window
               
                for j in range(Time_T*Time_window,Time_T*Time_window+Time_window+1):
                    temp_sum_1+=(expression_PPI_list[Dic_map[Node_PPI_1[Time_T][i]]][j]-temp_PPI_expression_avg_1)*(expression_PPI_list[Dic_map[Node_PPI_2[Time_T][i]]][j]-temp_PPI_expression_avg_2)
                    temp_sum_2+=(expression_PPI_list[Dic_map[Node_PPI_1[Time_T][i]]][j]-temp_PPI_expression_avg_1)**2
                    temp_sum_3+=(expression_PPI_list[Dic_map[Node_PPI_2[Time_T][i]]][j]-temp_PPI_expression_avg_2)**2                                 

                temp_r=temp_sum_1/((temp_sum_2*temp_sum_3)**0.5)
                temp_r=abs(temp_sum_1/((temp_sum_2*temp_sum_3)**0.5))
               
                if temp_r>1:
                    print("error!",temp_r)
                if temp_r>threshold_r:
                    
                    Adj_Matrix_Weight[Dic_map[Node_PPI_1[Time_T][i]],Dic_map[Node_PPI_2[Time_T][i]]]=temp_r
                    Adj_Matrix_Weight[Dic_map[Node_PPI_2[Time_T][i]],Dic_map[Node_PPI_1[Time_T][i]]]=temp_r
                else:
                    temp_num+=1
        df=pd.DataFrame(Adj_Matrix_Weight)
        df.to_excel("Adj_Matrix_Weight"+str(Time_T)+".xlsx")

  
