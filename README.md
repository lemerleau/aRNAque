# aRNAque (simple but efficient): A simple evolutionary program for and efficient RNA design.
<!--(@Author: [Nono Saha Cyrille Merleau](#) and [Matteo Smerlak](#) )-->

For a given target structure in a dot-bracket representation, the tool allows to generate good quality (low ED and MFE) RNA sequences with the corresponding strucure close to the input target. The method relies on local mutations of nucleotide and base pairs independently with respect to some probabilities: P_N and P_C.

![](images/mutation_example.png)

***Figure 1:** Mutation step illustration. (a) is a given target structure and (b) is a random compatible sequence from a population of RNA sequences. (c) is the mutated sequence where the non-base pair positions (in black color) are mutated independently of the base pair positions. One non-base pair position (3) and two base pair positions {(2,7);(11,16)} are mutated.*

The repo is organised as follows:
- [data](data/): The clean data used to produce the different plots presented in our paper.  The cleaned data are obtained by cleaning up the data generated from out file. for more details please refer to the python notebook [here](notebook/clean_data.ipynb)
- [docs] (docs): The files describing the evolutionary algorithm implemented.
- [images](images/): The plots (in pdf) used in the paper and the Python notebook code is in: [notebook/plots.ipynb](notebook/plots.ipynb).
- [src](src/): The source codes organized in three main parts:

    - [utility.py](src/utility.py): set of basic python functions usefull for our EA implementation.

    - [Landscape.py](src/Landscape.py): python class containing the information about the landscape to optimise.

    - [aRNAque.py](src/aRNAque.py): it contains the EA implementation, the initialization, mutation, selection and EA functions.

# Requirements
To be able to run aRNAque, the following softwares are required:

- [Python version 2.7](https://docs.anaconda.com/anaconda/user-guide/tasks/switch-environment/) or higher
- Numpy
- Pandas
- Scipy
- multiprocess: in case you would like to design many sequences in parallel.
- [ViennaRNA package](https://anaconda.org/bioconda/viennarna)
- python wrapper RNA (normally included in [ViennaRNA package](https://anaconda.org/bioconda/viennarna))
- seaborn: for the plots
- matplotlib: for the plotting part.

To install all requirements automatically, including the setup of a conda environment called `aRNAque` via miniconda, simply type the following command:

```
make requirements
```

The installation was tested on the following operating systems:

* MacOS Mojave
* Debian Xfce 4.12

# For the pseudoknotted RNA targets
## RNA folding programs for pseudoknotted targets:
  - IPknot ()
  - HotKnots ()
After installing the folding tools, make sure you have set an environmental variable for each tool HOTKNOTS_ROOT : for hotknots and IPKNOT for IPknot.

# How to run the program.
First, please clone the git repo using the command:

      $ git clone [repo link](#)
      $ cd aRNAque
      $ make requirements //In case the dependencies are not yet installed.  
      $ cd aRNAque/src/
      $ python aRNAque.py --target="((....)).((....)).((.....)).((....))"

For more details about the parameters please use:
```bash
❯ python aRNAque.py --help
usage: aRNAque.py [-h] [--target TARGET] [--job JOB] [-g G] [-n N] [-msf MSF]
                  [-sm SM] [-bp BP] [--Cs CS] [-EDg EDG] [-c C]
                  [--hairpin_boosting] [--folding_tool FOLDING_TOOL] [--log]
                  [--verbose] [--turner1999] [-seed SEED]

optional arguments:
  -h, --help            show this help message and exit
  --target TARGET, -t TARGET
                        Target RNA secondary structure in dot bracket
                        representation
  --job JOB, -j JOB     Number of EA runs (default: 1)
  -g G                  Number of generation (default: 150)
  -n N                  Population Size (default: 100)
  -msf MSF              maximum sequence found (default: 10)
  -sm SM                Selection method: the only possible values are {F,NED}
                        (default: NED)
  -bp BP                Distribution of nucleotide and base pairs. Possible
                        values are {GC,GC1,GC2,GC3,GC4,ALL}, please check the
                        online doc for more details (default: GC2)
  --Cs CS               sequence constraints: the lenght of the sequence
                        should be the same as the target. Example:
                        target=((....)), C=GNNNANNC (default: None)
  -EDg EDG              number of generation for Ensemble defect refinement
                        (default: 0)
  -c C                  Exponent of the zipf's distribution (default: None)
  --hairpin_boosting    Boost the hairpin loops. When false no hairpins
                        boosting (default: False)
  --folding_tool FOLDING_TOOL, -ft FOLDING_TOOL
                        folding tool to be used: v for RNAfold from viennarna,
                        l for LinearFold (default: v)
  --log                 Store the population for each instance of the inverse
                        folding in a folder (default: False)
  --verbose             Print the mean fitness evolution on a standard output
                        (default: False)
  --turner1999          Use the old energy parameters (default: False)
  -seed SEED            Seed for the initial population (default: None)
```
      $ python aRNAque.py --help

# For the pseudoknotted RNA targets
## thirdpart programs:
  - IPknot (for the folding of sequences)
  - HotKnots ()
After installing the folding tools make sure you have set an environmental variable for each tool.
