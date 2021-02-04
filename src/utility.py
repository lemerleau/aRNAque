
"""
    

"""
import pandas 
import random
import numpy 
import pp 
import os
import subprocess
import multiprocess
import time
import scipy
import argparse
import RNA
import uuid



def getHairepinCoord(sequence, target) : 
    cmd_input = str(sequence)+"\n"+str(target)+"\n"
    os.system("echo '"+cmd_input+"'| RNAeval -v |grep Hairpin | tr -d A-Z-a-z-'('-')' | cut -d ':' -f 1 > tmp/hairpin.loop")
    
    hairpins = []
    with open("tmp/hairpin.loop", "r") as loops : 
        data_r = loops.read().split("\n")
        for dt in data_r : 
            if dt != '' : 
                hairpins.append(tuple(numpy.array(list(map(int, dt.strip().split(','))))-1))
    loops.close()
    os.remove('tmp/hairpin.loop')

    return hairpins  


def boost_hairpins(sequence, coord) :
    stems = list(numpy.random.choice(["A"],len(range(coord[1]-coord[0]-1))))
    if coord[1]-coord[0]-1 >=4 : 
        stems[0] = "G"
        stems.insert(0,"G")
        stems.append("C")
        sequence[coord[0]:coord[1]+1] = stems
    else : 
        sequence[coord[0]+1:coord[1]] = stems
    return sequence


def pphamming(listOfStructures, landscape) : 
    pool = multiprocess.Pool(multiprocess.cpu_count())
    dists = pool.map(landscape.fitness,listOfStructures)
    pool.close()
    return dists

def ppens_defect(listOfSequences, landscape) : 
    pool = multiprocess.Pool(multiprocess.cpu_count())
    dists = pool.map(landscape.ens_defect,listOfSequences)
    pool.close()
    return dists

def get_bp_position(structure) : 
    position = RNA.ptable(structure)
    position = list(position[1:])

    base_paire_pos = []
    for i in range(len(position)) :
        if position[i] != 0 :
            if (position[i]-1,i) in base_paire_pos : 
                continue; 
            else : 
                base_paire_pos.append((i,position[i]-1))

    return base_paire_pos

#Logging population 
def bt_save_population(prev_pop, population,gen, root_path) : 
    data = []
    prev_data = []
    for i in range(len(population)) : 
        data.append([population[i].RNA_seq, population[i].RNA_structure,population[i].mfe, population[i].fitness])
        prev_data.append([prev_pop[i].RNA_seq, prev_pop[i].RNA_structure,prev_pop[i].mfe, prev_pop[i].fitness])

    
    dataFrame = pandas.DataFrame(data)
    prev_dataFrame = pandas.DataFrame(prev_data)
    prev_dataFrame.to_csv(root_path+"/prev_gen"+str(gen)+".csv")
    dataFrame.to_csv(root_path+"/gen"+str(gen)+".csv")


def save_population(population,gen, root_path, file_prefix) :
   
    population.to_csv(root_path+"/gen"+file_prefix+str(gen)+".csv")



#Evaluate the energy 
def ppeval(listOfSeqs, target, task) : 
    with open("tmp/rnaeval_in"+str(task), "w") as file_ : 
        for seq in listOfSeqs : 
            file_.write(seq.strip()+"\n"+target.strip()+"\n")
        file_.close()
          
    os.system("RNAeval -j --infile=tmp/rnaeval_in"+str(task)+" |tr -d A-Z,'(',')'|cut -d ' ' -f 2- > tmp/result_"+str(task))
    with open("tmp/result_"+str(task), "r") as file_ : 
        eval_ = file_.read().split()
    os.remove("tmp/result_"+str(task))
    os.remove("tmp/rnaeval_in"+str(task))
    return list(numpy.array(eval_, dtype=float))


def gen_point_mutation_dist(size,pos,c): 
    
    bp_pos = len(pos['bp_pos'])
    nbp_pos = len(numpy.where(numpy.array(pos["p_table"])==-1)[0])
    dist = {}
    if c !=None : 
        if 0<c<7.5:
            if bp_pos > 0 : 
                dist['bp'] = zipf_rvgen(1,bp_pos, bp_pos, size, c)
            else : 
                dist['bp'] = []
            dist['nbp'] = zipf_rvgen(1,nbp_pos, nbp_pos, size, c) 
        else : 
            dist['bp'] = modulo_dist(1,bp_pos, size) 
            dist['nbp'] = modulo_dist(1,nbp_pos,size)   
    else :
        if bp_pos > 0 : 
            dist['bp'] = numpy.ones(size, dtype=int)
        else : 
            dist['bp'] = []
        dist['nbp'] = numpy.ones(size, dtype=int)

    return dist


def ppfold(listOfSeqs,task) : 
    dataFrame= pandas.DataFrame(listOfSeqs)
    dataFrame.to_csv("tmp/sequences"+str(task),sep=" ",index=False, header=False)    
    os.system("RNAfold -j --infile=tmp/sequences"+str(task)+" --noPS > tmp/rnafold_result"+str(task)) #To change the energy par just add -P vrna185x.par  to RNAfold
    os.system("cut -d ' ' -f 1 tmp/rnafold_result"+str(task)+" | tr -d A-Z > tmp/result_"+str(task))
    os.system("cat tmp/rnafold_result"+str(task)+"|tr -d A-Z,'(',')' | cut -d ' ' -f 2- >tmp/mfes"+str(task))
    with open("tmp/result_"+str(task), "r") as file_ : 
        strc_ = file_.read().split()
        
    with open("tmp/mfes"+str(task), "r") as file_ : 
        mfes= file_.read().split()
     
    os.remove("tmp/result_"+str(task))
    os.remove("tmp/mfes"+str(task))
    os.remove("tmp/sequences"+str(task))
    os.remove("tmp/rnafold_result"+str(task))
 
    return strc_, list(numpy.array(mfes, dtype=float))
    