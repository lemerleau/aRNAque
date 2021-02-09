"""
    
    ######################################################################################################
                A SIMPLE EVOLUTIONARY ALGORITHM GUIDED BY LOCAL MUTATION FOR INVERSE RNA FOLDING
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
import Landscape 
import utility 


# Initialize a population of RNA.       
def init(pop_size, landscape, constraints=None) : 

    target = landscape.target 
    p_table = utility.get_bp_position(target)
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


#Mutation function 

def mutate(pop, pos, p_n, p_c,  distribution,constraints=None,loop_boosting=True) : 
    base_paire = ["GC","CG","AU","UA", "GU", "UG"]
    nucleotides = ["A", "G","U","C"]
    bp_pos = pos["bp_pos"]
    nbp_indexes = numpy.where(numpy.array(pos["p_table"])==-1)[0]
    mutants = []
    
    
    
    bp_points = distribution['bp']
    nb_points = distribution['nbp']
    
    for i in range(len(pop)) : 
        mutant = numpy.array(list(pop[i]))
        
	#Mutate a non base pair position
        nbp_indexes_to_mutate = numpy.random.choice(list(nbp_indexes), nb_points[i], replace=False)
        mutant[nbp_indexes_to_mutate]= numpy.random.choice(nucleotides, nb_points[i],p=p_n)


        #Mutate a base pair position
        if len(bp_pos) > 0 : 

            bp_choice = numpy.random.choice(base_paire,bp_points[i], p=p_c)
            bp_indexes_to_mutate = numpy.random.choice(range(len(bp_pos)), bp_points[i],replace=False)
            bp_to_mutate = numpy.array(bp_pos)[bp_indexes_to_mutate]
        
            for j in range(bp_points[i]) : 
                mutant[bp_to_mutate[j][0]] =bp_choice[j][0]
                mutant[bp_to_mutate[j][1]] =bp_choice[j][1]
            if loop_boosting : 
                for hp_pos in pos["hairpins"] : 
                    mutant = utility.boost_hairpins(mutant, hp_pos)

       
        if constraints!=None : 
            mutant[constraints['pos']] = constraints['sequence']
        
        mutants.append("".join(mutant))
        
    return mutants

# Selection function

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
def simple_EA(landscape, number_of_generation, init_pop, sm, log_folder,pos,p_n, p_c, msf=10, constraints=None,c=None, lb=True) : 
    """
    
    """
    print (" Starting of evolution ")
    prev_population = init_pop.copy(deep=True) #Initialize the population of RNA
    population_size =  len(init_pop)
    
    n = 0
    max_found = len(prev_population[prev_population["Fitness"]==str(1.0)])

    while (n<number_of_generation) and (max_found < msf):
           
        newgeneration = []
        best = elite(prev_population, 10, weight='Fitness')
        selected_ind = select(prev_population,population_size, sm)

        dist = utility.gen_point_mutation_dist(population_size,pos,c)      
        mutated = mutate(selected_ind,pos, p_n, p_c, dist,constraints,lb)+ list(best[:,0])

        newgeneration.append(mutated)

        structures, mfes = utility.ppfold(mutated,log_folder)
        newgeneration.append(structures)
        newgeneration.append(mfes)
        newgeneration.append(utility.pphamming(structures,landscape))
        evals = utility.ppeval(mutated, landscape.target,log_folder)
        newgeneration.append(evals)

        
        prev_population = pandas.DataFrame(numpy.array(newgeneration).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","Evals"])
        
        current_fitest = max(numpy.array(prev_population["Fitness"], dtype=float))

        n +=1
        max_found = len(prev_population[prev_population["Fitness"]==str(1.0)])
       
    return prev_population, n
 

def ED_refinement(landscape, number_of_generation,init_pop,log_folder,pos,p_n, p_c, constraints, c=None) : 
    
    print (" Starting ED refinement ")
    prev_population = init_pop.copy(deep=True) 

    
    population_size =  len(init_pop)
    n = number_of_generation
   
    while (n > 0) :

        newgeneration = []

        best = elite(prev_population, 10,weight='Fitness')
        max_fitness = max(numpy.array(prev_population['Fitness'], dtype=float)) 
        selected_ind = select(prev_population,population_size,'ED')
        
        dist = utility.gen_point_mutation_dist(population_size,pos,c) 
        mutated = mutate(selected_ind,pos, p_n, p_c,dist,constraints)+list(best[:,0])
        newgeneration.append(mutated)

        structures, mfes = utility.ppfold(mutated,log_folder)
        newgeneration.append(structures)
        newgeneration.append(mfes)
        newgeneration.append(utility.pphamming(structures,landscape))

        defects = utility.ppens_defect(mutated, landscape)
        newgeneration.append(defects)
        
        prev_population = pandas.DataFrame(numpy.array(newgeneration).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","ED"])
        seq_found = prev_population[prev_population["Fitness"]==str(max_fitness)]
        if len(seq_found.get_values()) > 0 :
            prev_population = seq_found
    

        n -=1
        
    return prev_population



def run(number_of_generation,pop, log_folder,landscape, p_n, p_c,msf,sm='NED',EDg=0, constraints=None, c=None,lb=True) : 
    target = landscape.target 
    p_table = utility.get_bp_position(target)
    
    evals = utility.ppeval(pop, target,log_folder)
    strcs, mfes = utility.ppfold(pop,log_folder)
    fitnesses = utility.pphamming(strcs, landscape)
    pos = {
        "p_table" : numpy.array(RNA.ptable(target)[1:]) - 1,
        "bp_pos"  : p_table,
        "hairpins": utility.getHairepinCoord(pop[0],target)
        }

    print (len(pos["hairpins"]), " Hairpin loop(s)")
    print (len(pos["bp_pos"])," loop(s) in total")
    init_pop = pandas.DataFrame(numpy.array([pop, strcs, mfes, fitnesses, evals]).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","Evals"])
    tic = time.time()
    best_pop, Mgen = simple_EA(landscape,number_of_generation,init_pop, sm,log_folder, pos, p_n, p_c, msf=msf, constraints=constraints,c=c)
    
    
    if EDg>0 : 
        ens_defects = utility.ppens_defect(list(best_pop["RNA_sequence"]), landscape)
        init_pop = pandas.DataFrame(numpy.array([list(best_pop["RNA_sequence"]), list(best_pop["RNA_structure"]), list(best_pop["Mfes"]), list(best_pop["Fitness"]), ens_defects]).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","ED"])
        
        best_pop = ED_refinement(landscape,EDg,init_pop,log_folder, pos, p_n, p_c,constraints,c)
        founds = best_pop[best_pop['Fitness']==str(max(numpy.array(best_pop['Fitness'], dtype=float)))]

        founds = founds.sort_values(by="ED", ascending=True)
        print (founds.get_values()[0])
    else : 
        founds = best_pop[best_pop['Fitness']==str(max(numpy.array(best_pop['Fitness'], dtype=float)))]
        founds = founds.sort_values(by='Mfes', ascending=False)
        print (founds.get_values()[0], Mgen)
    toc = time.time()
    print ("Time = ", toc-tic, "For method : ", sm)
    return (list(founds.get_values()[0]), Mgen)


def main() : 
    import sys
    
    if not os.path.exists("tmp/") : 
        os.mkdir('tmp')
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, argument_default=argparse.SUPPRESS)
    parser.add_argument('--target', type=str, help='Target RNA secondary structure in dot bracket representation')
    parser.add_argument('--job', type=int,default=1, help="Number of EA runs")
    parser.add_argument('-g', type=int, default=150, help="Number of generation")
    parser.add_argument('-n', type=int,default=100, help="Population Size")
    parser.add_argument('-msf', type=int,default=10, help="maximum sequence found")
    parser.add_argument('-sm', type=str, default='NED',help="Selection method: the only possible values are {F,NED}")
    parser.add_argument('-bp', type=str, default='GC2',help="Distribution of nucleotide and base pairs. Possible values are {GC,GC1,GC2,GC3,GC4,ALL}, please check the online doc for more details")
    parser.add_argument('--C', type=str, default=None,help="sequence constraints the lenght of the sequence should be the same as the target. Example: target=((....)), C=GNNNANNC")
    parser.add_argument('-EDg', type=int, default=0, help="number of generation for Ensemble defect refinement")
    parser.add_argument('-c', type=float,default=None, help="Exponent of the zipf's distribution" )
    parser.add_argument('-lb', type=bool,default=True, help="Boost the hairpin loops. When false no hairpins boosting" )
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
                print('Constraint sequence error: Please check the help using "python aRNAque.py --help"')
                
    else : 
        constraints = None
    landscape = Landscape.Landscape(target)
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
    job_server = pp.Server(30, ppservers=ppservers)
    pops = []
    pop = init(pop_size, landscape, constraints=None)
    
    for i in range(number_of_run) : 
        pops.append(pop)
     
    jobs = [(task , job_server.submit(run, (number_of_generation,pops[task], str(uuid.uuid4()),landscape, p_n, p_c,msf,sm,EDg,constraints,args.c,args.lb), (mutate,ED_refinement,elite, select, simple_EA),
     ("numpy", "RNA", "random","scipy.stats","pandas","os", "time","subprocess","multiprocess","utility","Landscape"))) for task in range(number_of_run)]

    result = []
    for task, job in jobs : 
        result.append(job())
        
    print (result)

    try : 
        os.removedirs('tmp')
    except OSError : 
        print("Error while cleaning the tmp directory files")
    
    


if __name__ == "__main__":
    main()
