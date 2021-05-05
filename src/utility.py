
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
import sys
sys.path.append('/usr/local/lib/python2.7/site-packages')
import RNA
import uuid
from folding_wrapper import *



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
    stems = ["A"]*(coord[1]-coord[0]-1)
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
    pk_pairs, bp, nbp= [] , [], []

    pairs = {
        'bp' : [],
        'pk' : [],
        'nbp' : []
    }

    for i,elt in enumerate(structure) :
        if elt =='[' :
            pk_pairs.append(i)
        elif elt ==']' :
            pairs['pk'] += [(pk_pairs.pop(),i)]

        elif elt == '(' :
            bp += [i]
        elif elt == ')' :
            pairs['bp'] += [(bp.pop(),i)]
        else :
            pairs['nbp'] += [i]
    return pairs

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
def ppeval(listOfSeqs, target) :
    task = uuid.uuid4()
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

def nthHarmonic(N,s) :

    harmonic = 1.00
    for i in range(2, N + 1) :
        harmonic += 1 / i**s

    return harmonic

def zipf_rvgen(low, high, N, size, a=6.) :
    choices = numpy.array(range(low, high+1))
    probs = choices**(-a) / nthHarmonic(N,a)
    return numpy.random.choice(choices,p=numpy.array(probs)/sum(probs),size=size)

def gen_point_mutation_dist(size,pos,c):

    bp_pos = len(pos['bp']+pos['pk'])
    nbp_pos = len(pos['nbp'])
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


def ppfold(listOfSeqs,tool) :

    if tool =="v" :
        return ppRNAfold(listOfSeqs)
    if tool =="l" :
        return pplinearFold(listOfSeqs)
    if tool=="c" :
        return ppcontextFold(listOfSeqs)
    if tool =="pk" :
        return ppKiss(listOfSeqs)
    if tool =="hk" :
        return ppHotknots(listOfSeqs)
