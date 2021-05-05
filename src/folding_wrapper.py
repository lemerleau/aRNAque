"""

"""

import os
import subprocess
import multiprocess as mp
import numpy
import sys
sys.path.append('/usr/local/lib/python2.7/site-packages')
import RNA
import uuid
import pandas
import array




def ppRNAfold(listOfSeqs) :
    task = uuid.uuid4()
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

def fold_with_contextfold(seq):

    p =subprocess.Popen(["java",'-cp',str(os.environ.get('CONTEXT_FOLD'))+"/bin","contextFold.app.Predict", "in:"+''.join(seq),"model:"+str(os.environ.get('CONTEXT_FOLD'))+"/trained/StHighCoHigh.model"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    rst, err = p.communicate()

    if err :
        print(err,"ERROR during the execution of contextFold.app.Predict program: please make sure the evironment variable CONTEXT_FOLD is set correctly")
        return;
    else :
        tmp_fc = RNA.fold_compound(seq)
        #print (str(rst).split())
        strc = str(rst).split()[1]
        mfe = tmp_fc.eval_structure(strc)
        assert len(strc)==len(seq)

        return (strc, mfe)


def ppcontextFold(listOfSeqs) :

    pool = mp.Pool(mp.cpu_count())
    result = numpy.array(pool.map(fold_with_contextfold, listOfSeqs))
    pool.close()

    return list(result[:,0]),list(result[:,1])


def fold_with_LinearFold(seq) :

    out_file = str(uuid.uuid4())
    cmd= "echo '"+seq+"' | "+os.environ.get("LINEAR_FOLD")+"/./linearfold > "+out_file
    os.system(cmd)
    with open(out_file) as file_:
        struct_dt = file_.read().split()

    file_.close()
    os.remove(out_file)
    return [struct_dt[1], float(struct_dt[-1][1:-1])]

def pKiss(seq) :

    cmd= "pKiss --strategy 'P' --mode='mfe' {} 2>/dev/null".format(seq)
    p = os.popen(cmd)
    rst = p.read().split()
    if len(rst) > 0 :
        return (rst[-1],rst[-2])
    else :
        print("ERROR during the folding with pKiss")
        return None

def hotknots(seq) :

    cmd= "./bin/HotKnots -s  {}".format(seq)
    p = os.popen(cmd)
    rst = p.read().split('\n')
    #print(rst)
    rst = rst[2].split()
    if len(rst) > 0 :
        #print (rst)
        return (rst[-2],rst[-1])
    else :
        print("ERROR during the folding with Hotknots")
        return None

def ppHotknots(listOfSeqs) :
    #pool = mp.Pool(mp.cpu_count())
    #result = numpy.array(pool.map(hotknots, listOfSeqs))
    #pool.close()
    result = []
    for s in listOfSeqs : 
        result.append(hotknots(s))
    result = numpy.array(result)
    return list(result[:,0]),list(result[:,1])


def ppKiss(listOfSeqs) :
    pool = mp.Pool(mp.cpu_count())
    result = numpy.array(pool.map(pKiss, listOfSeqs))
    pool.close()
    return list(result[:,0]),list(result[:,1])


def pplinearFold(listOfSeqs) :

    pool = mp.Pool(mp.cpu_count())
    result = numpy.array(pool.map(fold_with_LinearFold, listOfSeqs))
    pool.close()

    return list(result[:,0]),list(result[:,1])


def folding_with_mxfold(seq) :


    """
    def ppRNAfold(listOfSeqs) :
        task = uuid.uuid4()
        dataFrame= pandas.DataFrame(listOfSeqs)
        dataFrame.to_csv("tmp/sequences"+str(task),sep=" ",index=False, header=False)
        os.system("RNAfold -j --infile=tmp/sequences"+str(task)+" --noPS > tmp/rnafold_result"+str(task)) #To change the energy par just add -P vrna185x.par  to RNAfold

        with open("tmp/rnafold_result"+str(task)) as file_:
            fold_data = file_.read()
            #fold_data =np.array(fold_data).reshape(len(listOfSeqs),4)

        print ("".join(fold_data))
        #print (fold_data[:,-1])
        os.remove("tmp/rnafold_result"+str(task))

        return
    """



#Just for test
def main():

    #pseq = RNA.random_string(35,"AUGC")
    #print(pseq)
    #print(pKiss(pseq))


    if not os.path.exists("tmp/") :
        os.mkdir('tmp')
    random_seq = []

    for s in range(10):
        random_seq.append("".join(numpy.random.choice(["A","G","U", "C"], 30)))

    #rnastrand_sequences = pandas.read_csv("../data/clean_rnastrand_data.csv")["Sequence"].values

    #print (rnastrand_sequences)
    print ppHotknots(random_seq)
    #print (rst)
    #pandas.DataFrame(numpy.array([rnastrand_sequences,rst, mfes]).T, columns=['sequence','structure','mfe']).to_csv('../data/vrna_strand.csv')


if __name__ =="__main__" :
    main()
