"""
     @author : nonosaha@mis.mpg.de/cyrillecardinale@gmail.com

"""
#import libraries
import numpy as np
import RNA
import ast
import matplotlib.pyplot as plt
import scipy.stats as st
from pandas import read_csv
import os
import json


def entropy_seq(sequences) :
    p_k = []
    length = len(sequences[0])
    for n in ['A','C','G','U'] :
        p_n = []
        for i in range(length) :

            p = 0.
            for seq in sequences :
                if seq[i] == n :
                    p += 1.
            p_n.append(p/len(sequences))
        p_k += [p_n]
    print(st.entropy(p_k, axis=0))
    return sum(st.entropy(p_k, axis=0))


def main() :

    with open("../../data/diversity/distinct_stuctures.json") as file_ :
        distinct_strucs = json.load(file_)
        file_.close()

    with open("../../data/diversity/distinct_sequences.json") as file_ :
        distinct_seq = json.load(file_)
        file_.close()

    with open("../../data/diversity/entropies.json") as file_ :
        diversities = json.load(file_)
        file_.close()

    with open("../../data/diversity/mean_fitness.json") as file_ :
        mean_fitness = json.load(file_)
        file_.close()

    with open("../../data/diversity/max_fitness.json") as file_ :
        max_fitness = json.load(file_)
        file_.close()


    figure = plt.figure(constrained_layout=True, figsize=(9,6))
    gs = figure.add_gridspec(nrows=2, ncols=2, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    plt.title("(A)", fontsize=15)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.ylabel(r"Mean entropy ($<H_t>$)", fontsize=12)
    plt.xlabel(r'Generation time ($t$)',fontsize=12)

    for f in range(10):
        if f == 0 :
            plt.plot(diversities['op'][str(f)] , color="darkorange", label="One point mutation")
            plt.plot(diversities['levy'][str(f)] , color="deepskyblue", label="Lévy mutation")
        plt.plot(diversities['op'][str(f)] , color="darkorange")
        plt.plot(diversities['levy'][str(f)] , color="deepskyblue")
    plt.legend()

    ax = figure.add_subplot(gs[0,1])
    plt.title("(B)", fontsize=15)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    for f in range(10):
        if f == 0 :
            plt.plot(max_fitness['op'][str(f)],diversities['op'][str(f)] ,color="darkorange",label='One point mutation')
            plt.plot( max_fitness['levy'][str(f)],diversities['levy'][str(f)] , color='deepskyblue', label='Levy mutation')

        plt.plot(max_fitness['op'][str(f)],diversities['op'][str(f)] ,color="darkorange")
        plt.plot( max_fitness['levy'][str(f)],diversities['levy'][str(f)], color='deepskyblue')

    #ax.set(yscale='log')
    #ax.set(xscale='log')
    ax.set_ylabel(r"Mean entropy ($<H_t>$)", fontsize=12)
    ax.set_xlabel(r"Max fitness ($F_t$)", fontsize=12)
    #plt.legend()



    ax = figure.add_subplot(gs[1,0])

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.title("(C)", fontsize=15)
    for f in range(10) :

        if f == 0 :
            plt.plot(max_fitness['op'][str(f)], color='darkorange', label='One point mutation')
            plt.plot(max_fitness['levy'][str(f)], color='deepskyblue', label='Levy mutation')

        plt.plot(max_fitness['op'][str(f)], color='darkorange')
        plt.plot(max_fitness['levy'][str(f)], color='deepskyblue')


    #plt.legend()
    plt.ylabel("Max fitness", fontsize=12 )
    plt.xlabel(r'Generation time($t$)',fontsize=12)

    ax2 = plt.axes([0.25, 0.13, 0.2, 0.2])
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    for f in range(10) :
        print("folder ", f)

        if f == 0 :
            plt.plot(mean_fitness['op'][str(f)], color='darkorange', label='One point mutation')
            plt.plot(mean_fitness['levy'][str(f)], color='deepskyblue', label='Levy mutation')

        plt.plot(mean_fitness['op'][str(f)], color='darkorange')
        plt.plot(mean_fitness['levy'][str(f)], color='deepskyblue')


    ax = figure.add_subplot(gs[1,1])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.title("(D)", fontsize=15)
    for f in range(10) :

        if f== 1 :
            plt.plot(distinct_seq["op"][str(f)],distinct_strucs['op'][str(f)], color='darkorange', label="One point mutation")
            plt.plot(distinct_seq["levy"][str(f)],distinct_strucs['levy'][str(f)], color='deepskyblue', label="Levy mutation")

        plt.plot(distinct_seq["op"][str(f)],distinct_strucs['op'][str(f)], color='darkorange')
        plt.plot(distinct_seq["levy"][str(f)],distinct_strucs['levy'][str(f)], color='deepskyblue')


    plt.xlabel(r"Distinct sequences",fontsize=12)
    plt.ylabel(r"Distinct structures", fontsize=12)
    plt.savefig("../../images/diversity.pdf")
    plt.show()



if __name__ == '__main__' :
    main()
