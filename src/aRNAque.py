"""
    @author: Nono Saha Cyrille Merleau  
    @email: nonosaha@mis.mpg.de


    ######################################################################################################
                A SIMPLE EVOLUTIONARY ALGORITHM INSPIRED BY LEVY-FLIGHT FOR INVERSE RNA FOLDING
    ######################################################################################################
"""

#Importting necessary python libraries
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

class Landscape (object) : 

    def __init__ (self, target) : 
        self.target = target 
    
    def fitness (self, structure) : 
        return 1./(1+RNA.hamming_distance(self.target, structure))

    def ens_defect(self, sequence) : 
        p = subprocess.Popen("defect", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cmd_input = sequence+"\n"+self.target+"\n"
        
        defect, error = p.communicate(cmd_input)
        defect = defect.split("\n")
        
        return 1/float(defect[-3])

    def ens_defect2(self, sequence) : 

        fc = RNA.fold_compound(sequence)
        fc.pf()
        fc.bpp()
        ed = fc.ensemble_defect(self.target)
        return ed

def getHairepinCoord(sequence, target) : 
    cmd_input = str(sequence)+"\n"+str(target)+"\n"
    os.system("echo '"+cmd_input+"'| RNAeval -v |grep Hairpin | tr -d A-Z-a-z-'('-')' | cut -d ':' -f 1 > tmp/hairpin.loop")
    
    hairpins = []
    with open("tmp/hairpin.loop", "r") as loops : 
        data_r = loops.read().split("\n")
        for dt in data_r : 
            if dt != '' : 
                hairpins.append(tuple(numpy.array(list(map(int, dt.strip().split(','))))-1))
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
    dists = pool.map(landscape.ens_defect2,listOfSequences)
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

#Mutation function 
def mutateOne(seq, mut_probs,target,pos, p_n, p_c, mut_bp=0.1) :  # best = 0.5, 0.1, 1/6
    base_paire = ["GC","CG","AU","UA", "GU", "UG"]
    nucleotides = ["A", "G","U","C"]
    p_table = pos["bp_pos"]
    RNA_seq = numpy.array(list(seq))
    r = numpy.random.rand(len(seq))
    mut_pos =RNA_seq[r<mut_probs] 
    choices = numpy.random.choice(nucleotides, len(mut_pos), p=p_n)
    RNA_seq[r<mut_probs] = choices 
    apply = []
    for bp_cord in p_table : 
        r = random.uniform(0,1)
        
        if r < mut_bp : 
            bp = numpy.random.choice(base_paire,1, p=p_c)
            
            RNA_seq[bp_cord[0]] = bp[0][0]
            RNA_seq[bp_cord[1]] = bp[0][1]
        
        if bp_cord in pos["hairpins"] : 
                RNA_seq = boost_hairpins(RNA_seq,target, bp_cord)

    return ''.join(RNA_seq)


def mutateAll(population, mut_prob, target,pos, p_n, p_c) : 
    mutated_pop = [] 
    for individual in population : 
        mutated_pop.append(mutateOne(individual,mut_prob,target,pos,p_n, p_c))
    return mutated_pop


def mutateOnePoint(seq): 
    nucleotides = ["A", "C", "G" , "U"]
    mutate_pos = numpy.random.randint(low=0, high=len(seq))
    RNA_seq = list(seq)
    RNA_seq[mutate_pos] = numpy.random.choice(nucleotides, 1)[0]
    return "".join(RNA_seq)


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

def modulo_dist(low, high, size) : 
    op = numpy.ceil(size*0.8)
    hpm= size-op 
    
    A = numpy.ones(int(op), dtype=int)
    B = numpy.random.randint(int(high/2.),high, size=int(hpm)) 
    
    return numpy.concatenate((A,B))

def mutate(pop, pos, p_n, p_c,  distribution,constraints=None) : 
    base_paire = ["GC","CG","AU","UA", "GU", "UG"]
    nucleotides = ["A", "G","U","C"]
    bp_pos = pos["bp_pos"]
    nbp_indexes = numpy.where(numpy.array(pos["p_table"])==-1)[0]
    mutants = []
    
    
    
    bp_points = distribution['bp']
    nb_points = distribution['nbp']
    
    for i in range(len(pop)) : 
        mutant = numpy.array(list(pop[i]))
        
        #Mutate a base pair position
        if len(bp_pos) > 0 : 

            bp_choice = numpy.random.choice(base_paire,bp_points[i], p=p_c)
            bp_indexes_to_mutate = numpy.random.choice(range(len(bp_pos)), bp_points[i],replace=False)
            bp_to_mutate = numpy.array(bp_pos)[bp_indexes_to_mutate]
        
            for j in range(bp_points[i]) : 
                mutant[bp_to_mutate[j][0]] =bp_choice[j][0]
                mutant[bp_to_mutate[j][1]] =bp_choice[j][1]
       #Mutate a non base pair position
        nbp_indexes_to_mutate = numpy.random.choice(list(nbp_indexes), nb_points[i], replace=False)
        mutant[nbp_indexes_to_mutate]= numpy.random.choice(nucleotides, nb_points[i],p=p_n)
        if constraints!=None : 
            mutant[constraints['pos']] = constraints['sequence']
        
        mutants.append("".join(mutant))
        
    return mutants

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
    

def select(population, size, method='F') : 
    
    if method == 'F' : 
        weights  = numpy.array(population["Fitness"], dtype = float)
    elif method == 'NED' : 
        evals = numpy.array(population["Evals"], dtype = float)
        mfes  = numpy.array(population["Mfes"], dtype = float)
    
        delta_mfe = evals - mfes 
        delta_mfe = delta_mfe / sum(delta_mfe)
        weights = (1-delta_mfe)**100
    elif method =='ED' : 
        weights = 1./numpy.array(population["ED"], dtype = float)
    else : 
        print ('The method '+method+ 'does not exist')
        return;
    
    selected = numpy.random.choice(list(population["RNA_sequence"]),size=size,p=weights/sum(weights))
    return list(selected)


#Elitism method that aims to copy the number of fittest sequences to the next generation 
def elite(population, size, weight='Fitness') : 
    sorted_pop = population.sort_values(by=weight, ascending=False)
                
    return sorted_pop.get_values()[:size]
 


#Main Evolutionary algorithm
def simple_EA(landscape, number_of_generation, mut_probs, init_pop, sm, log_folder,pos,p_n, p_c, msf=10, constraints=None,c=None) : 
    """
    if c !=None : 
        if c in [3.5, 5.5] : 
            root_folder = "logs/mutation/zipf/"+str(c)+"/"+log_folder+"/"
        else : 
            root_folder ="logs/mutation/modal/"+log_folder+"/"
    else : 
        root_folder = "logs/mutation/onepoint/"+log_folder+"/"
    #os.makedirs(root_folder) 
    """
    print (" Starting of evolution ")
    prev_population = init_pop.copy(deep=True) #Initialize the population of RNA
    population_size =  len(init_pop)
    
    n = 0
    #save_population(prev_population,n,root_folder,"")
    max_found = len(prev_population[prev_population["Fitness"]==str(1.0)])

    while (n<number_of_generation) and (max_found < msf):
           
        newgeneration = []
        best = elite(prev_population, 10, weight='Fitness')
        selected_ind = select(prev_population,population_size, sm)

        dist = gen_point_mutation_dist(population_size,pos,c)      
        mutated = mutate(selected_ind,pos, p_n, p_c, dist,constraints)+ list(best[:,0])

        newgeneration.append(mutated)

        structures, mfes = ppfold(mutated,log_folder)
        newgeneration.append(structures)
        newgeneration.append(mfes)
        newgeneration.append(pphamming(structures,landscape))
        evals = ppeval(mutated, landscape.target,log_folder)
        newgeneration.append(evals)

        
        prev_population = pandas.DataFrame(numpy.array(newgeneration).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","Evals"])
        
        current_fitest = max(numpy.array(prev_population["Fitness"], dtype=float))

        n +=1
        max_found = len(prev_population[prev_population["Fitness"]==str(1.0)])
        #save_population(prev_population,n,root_folder,"")
       
    return prev_population, n
 

def ED_refinement(landscape, number_of_generation,init_pop,log_folder,pos,p_n, p_c, constraints, c=None) : 
    
    print (" Starting ED refinement ")
    prev_population = init_pop.copy(deep=True) #Initialize the population of RNA

    
    population_size =  len(init_pop)
    n = number_of_generation
   
    while (n > 0) :

        newgeneration = []

        best = elite(prev_population, 10,weight='Fitness')
        max_fitness = max(numpy.array(prev_population['Fitness'], dtype=float)) 
        selected_ind = select(prev_population,population_size,'ED')
        
        dist = gen_point_mutation_dist(population_size,pos,c) 
        mutated = mutate(selected_ind,pos, p_n, p_c,dist,constraints)+list(best[:,0])
        newgeneration.append(mutated)

        structures, mfes = ppfold(mutated,log_folder)
        newgeneration.append(structures)
        newgeneration.append(mfes)
        newgeneration.append(pphamming(structures,landscape))

        defects = ppens_defect(mutated, landscape)
        newgeneration.append(defects)
        
        prev_population = pandas.DataFrame(numpy.array(newgeneration).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","ED"])
        seq_found = prev_population[prev_population["Fitness"]==str(max_fitness)]
        if len(seq_found.get_values()) > 0 :
            prev_population = seq_found
    

        n -=1
        
    return prev_population
       
def init(pop_size, landscape, constraints=None) : 

    target = landscape.target 
    p_table = get_bp_position(target)
    init_depth =len(target)
   
    nucleotides = ["A", "U" , "G", "C"]
    base_paire = ["GC","CG","AU","UA", "GU", "UG"]
    pop = []
    i = 0
    
    while i < pop_size : 
        
        if i < 4: 
            arn = numpy.random.choice(nucleotides[i:i+1],init_depth)
        else : 
            arn = numpy.random.choice(nucleotides,len(target))
            
        for bp_cord in p_table : 
            bp = numpy.random.choice(base_paire,1)
            arn[bp_cord[0]] = bp[0][0]
            arn[bp_cord[1]] = bp[0][1]
        if constraints!=None : 
                arn[constraints['pos']] = constraints['sequence']
        pop.append(''.join(arn))
        i = len(pop)

    return pop

def run(number_of_generation,pop, mut_probs, log_folder,landscape, p_n, p_c,msf,sm='NED',EDg=0, constraints=None, c=None) : 
    target = landscape.target 
    p_table = get_bp_position(target)
    
    evals = ppeval(pop, target,log_folder)
    strcs, mfes = ppfold(pop,log_folder)
    fitnesses = pphamming(strcs, landscape)
    pos = {
        "p_table" : numpy.array(RNA.ptable(target)[1:]) - 1,
        "bp_pos"  : p_table,
        "hairpins": getHairepinCoord(pop[0],target)
        }

    print (len(pos["hairpins"]), " Hairpin loop(s)")
    print (len(pos["bp_pos"])," loop(s) in total")
    init_pop = pandas.DataFrame(numpy.array([pop, strcs, mfes, fitnesses, evals]).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","Evals"])
    tic = time.time()
    best_pop, Mgen = simple_EA(landscape,number_of_generation, mut_probs,init_pop, sm,log_folder, pos, p_n, p_c, msf=msf, constraints=constraints,c=c)
    
    
    if EDg>0 : 
        ens_defects = ppens_defect(list(best_pop["RNA_sequence"]), landscape)
        init_pop = pandas.DataFrame(numpy.array([list(best_pop["RNA_sequence"]), list(best_pop["RNA_structure"]), list(best_pop["Mfes"]), list(best_pop["Fitness"]), ens_defects]).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","ED"])
        
        best_pop = ED_refinement(landscape,EDg,init_pop,log_folder, pos, p_n, p_c,constraints,c)
        founds = best_pop[best_pop['Fitness']==str(max(numpy.array(best_pop['Fitness'], dtype=float)))]

        founds = founds.sort_values(by="ED", ascending=True)
        #print (founds.get_values()[0])
    else : 
        founds = best_pop[best_pop['Fitness']==str(max(numpy.array(best_pop['Fitness'], dtype=float)))]
        founds = founds.sort_values(by='Mfes', ascending=False)
        #print (founds.get_values()[0], Mgen)
    toc = time.time()
    print ("Time = ", toc-tic, "For method : ", sm)
    return (list(founds.get_values()[0]), Mgen)


def main() : 
    import sys
    
    if not os.path.exists("tmp/") : 
        os.mkdir('tmp')
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, argument_default=argparse.SUPPRESS)
    parser.add_argument('--target', type=str, help='Target RNA secondary structure in dot bracket')
    parser.add_argument('--job', type=int,default=1, help="Number of runs")
    parser.add_argument('-g', type=int, default=150, help="Number of generation")
    parser.add_argument('-n', type=int,default=100, help="Population Size")
    parser.add_argument('-msf', type=int,default=10, help="maximum sequence found")
    parser.add_argument('-sm', type=str, default='NED',help="Selection method: the only possible values are {F,NED}")
    parser.add_argument('-bp', type=str, default='GC',help="base pair constraints possible values are {GC,GC1,ALL}")
    parser.add_argument('--C', type=str, default=None,help="sequence constraints the lenght of the sequence should be the same as the target. Example: target=((....)), C=GNNNANNC")
    parser.add_argument('-EDg', type=int, default=0, help="number of generation for Ensemble defect refinement")
    parser.add_argument('-c', type=float,default=None, help="Exponent of the zipf's distribution" )
    args = parser.parse_args()
    target = args.target
    sm = args.sm
    msf = args.msf
    bp = args.bp 
    main_sequence = args.C
    EDg = args.EDg 
    constraints = {}
    if main_sequence !=None : 
        if len(main_sequence) == len(target) : 
            if set(list(main_sequence)).issubset(set(["A","C","G","U","N"])) : 
                constraints["pos"] = numpy.where(numpy.array(list(main_sequence))!="N")
                constraints["sequence"] = numpy.array(list(main_sequence))[constraints["pos"]]
                print(constraints)
            else : 
                print('Constraint sequence error: Please check the help using "python rnaevol.py --help"')
                
    else : 
        constraints = None
    landscape = Landscape(target)
    print ("=====================================================================================================")
    print ("Solving for the target = ", target )
    print ("=====================================================================================================")
    number_of_run = args.job
    init_depth =len(target)
    mut_prob = 1./init_depth
    number_of_generation = args.g
    pop_size = args.n
     
    if bp == 'GC5' : 
        p_n = [0.25,0.65,0.05,0.05] 
        p_c = [0.4,0.5,0.1,0.0,0.0,0.0] 
    if bp == 'GC' : 
        p_n = [0.25,0.25,0.25,0.25] 
        p_c = [0.5,0.5,0.0,0.0,0.0,0.0] 
    if bp == 'GC1' : 
        p_n = [0.7,0.1,0.1,.1] 
        p_c = [0.3, 0.2, 0.2, 0.1,0.1,0.1] 
    if bp == 'GC4' : 
        p_n = [0.95,0.0,0.05,0.0] 
        p_c = [0.4, 0.4, 0.2, 0.,0.,0.] 
    if bp == 'GC3' : 
        p_n = [0.75,0.1,0.1,.05] 
        p_c = [0.4, 0.5, 0.1, 0.,0.,0.] 
    if bp == 'GC2' : 
        p_n = [0.7,0.1,0.1,.1] 
        p_c = [0.4, 0.5, 0.1, 0.,0.,0.] 
    if bp == 'ALL' : 
        p_n = [0.25,0.25,0.25,.25] #default = [0.25,0.25,0.25,.25] [0.25,0.65,0.05,.05] [0.7,0.1,0.1,.1] ["A", "G","U","C"]
        p_c = [0.2,0.2,0.1,0.1,0.2,0.2]  #[0.2,0.2,0.1,0.1,0.2,0.2] #[0.4, 0.5, 0.1, 0.,0.,0.] ["GC","CG","AU","UA", "GU", "UG"]
    ppservers = ()
    mut_probs = numpy.array(RNA.ptable(target)[1:])
    mut_probs = mut_probs + mut_prob
    mut_probs[mut_probs>mut_prob] = 0

    job_server = pp.Server(30, ppservers=ppservers)
    pops = []
    pop = init(pop_size, landscape, constraints=None)
    
    for i in range(number_of_run) : 
        pops.append(pop)
     
    
    print("============ Zipf point mutation ==============")
    jobs = [(task , job_server.submit(run, (number_of_generation,pops[task], mut_probs, str(uuid.uuid4()),landscape, p_n, p_c,msf,sm,EDg,constraints,args.c), (mutate,gen_point_mutation_dist,ED_refinement,elite, select,ppeval,ppfold, pphamming,mutateAll,get_bp_position,getHairepinCoord, boost_hairpins, mutateOne, modulo_dist, save_population, zipf_rvgen,nthHarmonic, simple_EA, ppens_defect),("numpy", "RNA", "random","scipy.stats","pandas","os", "time","subprocess","multiprocess"))) for task in range(number_of_run)]

    result = []
    for task, job in jobs : 
        result.append(job())
        
    print (result)

if __name__ == "__main__":
    main()
