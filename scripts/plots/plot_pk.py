import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb






def main() :
    #Load pk data
    ipknot_df = pd.read_csv("../../data/PseudoBase++/result_ipknot.csv")
    new_data = []
    s = 20
    labels= { "["+str(t)+ "-" +str(t+s)+"]" : (t, t+s, []) for t in range(min(ipknot_df['Length'].values),max(ipknot_df['Length'].values),s)}

    print(labels)

    for key in labels.keys() :
        for l in range(labels[key][0],labels[key][1]):
            val = ipknot_df[ipknot_df["Length"]==l].values.tolist()
            if len(val) > 0 :
                for elt in val :
                    nl = list(elt)
                    nl [-1] = nl[-1] + 1
                    new_data += [nl+[key]]
    print(len(ipknot_df), len(new_data))
    new_dt = pd.DataFrame(new_data, columns=["id","sequence","structure","target", "Length","Hamming distance",'bp_density'
    ,'Mutation mode','Number of generations', "Length group"])
    op_df = ipknot_df[ipknot_df["Mutation mode"]=='OP']["Generation"].values.reshape(253,20)
    op_med = [np.median(df) for df in op_df]

    levy_df = ipknot_df[ipknot_df["Mutation mode"]=='Levy']["Generation"].values.reshape(253,20)
    levy_med = [np.median(df) for df in levy_df]


    figure = plt.figure(constrained_layout=True, figsize=(9,6))
    gs = figure.add_gridspec(nrows=2, ncols=1, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    #sb.set(rc={'figure.figsize': (20., 8.27)})
    #plt.ylabel('Hamming Distance')
    plt.title("(A)")
    plt.xticks(rotation=90, fontsize=12)
    plt.title("(A)", fontsize=15)
    ax.set_ylabel('Hamming distance',fontsize=12)
    #ax.set_xlabel('Length Group', fontsize=12)
    sb_bx = sb.boxplot(ax=ax, y='Hamming distance', x='Length group', hue='Mutation mode', data=new_dt, palette={
    "Levy": "deepskyblue",
    "OP": "darkorange"
    })
    sb_bx.set(xticklabels=[])
    sb_bx.set(xlabel=None)


    ax = figure.add_subplot(gs[1,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.xticks(rotation=90)
    plt.title("(B)", fontsize=15)
    ax.set_ylabel('Hamming distance', fontsize=12)
    ax.set_xlabel('Length group', fontsize=12)
    ax.set(yscale="log")
    sb.boxplot(ax=ax, y='Number of generations', x='Length group', hue='Mutation mode', data=new_dt,palette={
    "Levy": "deepskyblue",
    "OP": "darkorange"
    })
    plt.legend([],[], frameon=False)
    plt.savefig('../../images/PseudoBase++/pkbase_ipknotV2.pdf')
    plt.show()

    diff = [op_med[i]-levy_med[i] for i in range(253)]

    figure = plt.figure(constrained_layout=True, figsize=(7,5))
    gs = figure.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.title("(A)",fontsize=15)
    plt.xlabel(r"Median number of generations ($t$)")
    plt.hist(levy_med, bins=10, label="Levy", color="deepskyblue")
    plt.hist(op_med, bins=10, alpha=0.3,label="OP", color="darkorange")
    #plt.legend()

    op_success = [(20-list(gens).count(200))/0.2 for gens in op_df]
    levy_success = [(20-list(gens).count(200))/0.2 for gens in levy_df]
    ax = figure.add_subplot(gs[0,1])

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.title("(B)",fontsize=15)
    plt.xlabel("Success rate (%)")
    plt.hist(levy_success, bins=10, label="Levy", color="deepskyblue")
    plt.hist(op_success, bins=10, alpha=0.4,label="OP", color="darkorange")
    plt.legend()
    plt.savefig("../../images/PseudoBase++/pk_histo.pdf")
    plt.show()

    print()



if __name__=="__main__" :
    main()