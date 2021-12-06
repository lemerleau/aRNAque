from matplotlib import pyplot as plt
import numpy as np

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS



def nthHarmonic(N,s) :

    harmonic = 1.00
    for i in range(2, N + 1) :
        harmonic += 1 / i**s

    return harmonic

def zipf_rvgen(low, high, N, size, a=6.) :
    choices = np.array(range(low, high+1))
    probs = choices**(-a) / nthHarmonic(N,a)
    return np.random.choice(choices,p=np.array(probs)/sum(probs),size=size)

def bino_gen(mu, length, size) :
    rst = [list(np.random.binomial(1, mu, length)).count(1) for i in range(size)]
    return rst

def parseArguments():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, argument_default=SUPPRESS)
    parser.add_argument("--length", '-l' ,type=int, help= "Lenght of the Sequence")
    parser.add_argument('-c', type=float, default=1.5, help='Levy of zipf law exponent')
    parser.add_argument('-ps', type=int, default=100, help='Population size')
    return parser.parse_args()

def main() :
    args = parseArguments()
    zipf_dist = list(zipf_rvgen(1,args.length, args.length, args.ps, args.c))

    barplot_data = []

    for i in set(zipf_dist) :
        barplot_data +=[[i,zipf_dist.count(i)]]


    barplot_data = np.array(barplot_data)
    p = barplot_data[:,1]/sum(barplot_data[:,1])
    mean_zipf = sum(p*barplot_data[:,0])
    print(f"The mean of the Zipf's distribution is {mean_zipf}")
    mu = mean_zipf/args.length
    print(f"The Binomial mutation rate is {mu}")


    bino_dist = bino_gen(mu, args.length, args.ps)
    barplot_data_bino = []

    for i in set(bino_dist) :
        barplot_data_bino +=[[i,bino_dist.count(i)]]

    barplot_data_bino = np.array(barplot_data_bino)
    figure = plt.figure(constrained_layout=True, figsize=(6,4))
    gs = figure.add_gridspec(nrows=1, ncols=1, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    #plt.title("Binomial and Zipf distribution centered at the mean = "+str(np.round(mean_zipf,2)))
    plt.xlabel("Number of point mutation", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.bar(barplot_data_bino[:,0], barplot_data_bino[:,1], width=.8, align='center', log=True, label=r"Binomial with $\mu="+str(np.round(mu,3))+"$", color="darkorange")
    plt.bar(barplot_data[:,0], barplot_data[:,1], width=.8, align='center', log=True, label=r"Zipf with $c="+str(args.c)+"$", alpha=0.5, color="deepskyblue")
    plt.legend()
    plt.savefig("../../images/bino_zipf.pdf")
    plt.show()




if __name__ =="__main__" :
    main()
