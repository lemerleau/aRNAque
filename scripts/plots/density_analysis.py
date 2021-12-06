"""
     @author : nonosaha@mis.mpg.de/cyrillecardinale@gmail.com

"""
#import libraries
import numpy as np
import RNA
import ast
import matplotlib.pyplot as plt
import scipy.stats as st
import pandas as pd
import os
import json
import seaborn as sb


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


    with open("../../data/Eterna100/V1/eterna_ds.json") as js_file :
        eterna_density = json.load(js_file)
        js_file.close()

    eterna_ds = {
        'low' : [],
        'high': [],
    }
    eterna_ds_Y = [0., 0.]
    for key in eterna_density.keys() :
        if float(key) <= 0.5 :
            eterna_ds_Y[0] += len(eterna_density[key])
            eterna_ds['low'] += eterna_density[key]
        else :
            eterna_ds_Y[1] += len(eterna_density[key])
            eterna_ds['high'] += eterna_density[key]

    pk_ds = {}

    pk_df = pd.read_csv("../../data/PseudoBase++/pkbase.csv")
    for tg in pk_df.values[:,-1] :
        pk_ds[tg] = np.round((tg.count('(')*2. + tg.count('[')*2. + tg.count('{')*2.)/len(tg),1)
    pk_ds_X = set(list(pk_ds.values()))

    pkbp_ds = list(pk_ds.values())

    pk_ds_Y = [0., 0.]

    for key in pk_ds.keys() :

        if pk_ds[key] <= 0.5 :
            pk_ds_Y[0] += 1
        else :
            pk_ds_Y[-1] += 1

    eterna_1999GC2 = pd.read_csv('../../data/Eterna100/V1/eterna1999_zipf_GC2.csv').values[:,-1].reshape(100,5)
    targets = []

    with open("../../data/Eterna100/V1/eterna_1999_zipfGC2.log") as file_ :
        lines = file_.readlines()

        for line in lines :
            if line.startswith("('Solving") :
                targets += [ast.literal_eval(line)[-1]]

    eterna199_levy_ds = {
        'Low' : [],
        'High' : []
    }

    for i,t in enumerate(targets) :
        if t in eterna_ds['low'] :
            eterna199_levy_ds["Low"] += eterna_1999GC2[i].tolist()
        else:

            eterna199_levy_ds["High"] += eterna_1999GC2[i].tolist()

    plot_df = []
    for key in eterna199_levy_ds.keys() :
        for val in eterna199_levy_ds[key] :
            plot_df += [[key, val+1,"Levy"]]

    eterna_1999GC2_OP = pd.read_csv('../../data/Eterna100/V1/eterna1999op_GC2.csv').values[:,-1].reshape(100,5)

    eterna199_OP_ds = {
    'Low' : [],
    'High' : []
    }

    for i,t in enumerate(targets) :
        if t in eterna_ds['low'] :
            eterna199_OP_ds["Low"] += eterna_1999GC2_OP[i].tolist()
        else:

            eterna199_OP_ds["High"] += eterna_1999GC2_OP[i].tolist()

    for key in eterna199_OP_ds.keys() :
        for val in eterna199_OP_ds[key] :
            plot_df += [[key, val+1,"OP"]]
    df_plot = pd.DataFrame(plot_df, columns=['Base Pair density', r"Generation ($t$)", "Mutation type"])

    
    figure = plt.figure(constrained_layout=True, figsize=(10,4))
    gs = figure.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax.set_xlabel("Base Pair density", fontsize=12)
    ax.set_ylabel(r"Generation ($t$)", fontsize=12)
    ax.set(yscale="log")
    plt.title("(A)", fontsize=15)
    sb.boxplot(data=df_plot, x="Base Pair density", y=r"Generation ($t$)", hue="Mutation type",ax=ax, palette={
        "Levy": "deepskyblue",
        "OP": "darkorange"
    })
    plt.legend(loc="upper left")
    ax = figure.add_subplot(gs[0,1])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.title("(B)", fontsize=15)
    plt.xlabel('Base pair density', fontsize=12)
    plt.ylabel('# of targets', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.bar(["Low","High"], height=pk_ds_Y, align='center', width=0.3, label='Pseudobase++',color='wheat')
    for i in range(2):
        plt.annotate(str(np.round(pk_ds_Y[i]*100/sum(pk_ds_Y), 2))+"%", xy=(i,pk_ds_Y[i]-((pk_ds_Y[i]*100/sum(pk_ds_Y))-13)), ha='center', va='bottom')
    plt.bar(["Low","High"], height=eterna_ds_Y, align='center', width=0.3, label='Eterna100', color='tan')
    for i in range(2):
        plt.annotate(str(np.round(eterna_ds_Y[i]*100/sum(eterna_ds_Y), 2))+"%", xy=(i,eterna_ds_Y[i]-35), ha='center', va='bottom')

    plt.legend(loc="upper left")
    plt.savefig('../../images/levy_analysis.pdf')
    plt.show()

if __name__ == '__main__' :
    main()
