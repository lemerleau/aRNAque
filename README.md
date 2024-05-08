# aRNAque.
<!--(@Author: [Nono Saha Cyrille Merleau](#) and [Matteo Smerlak](#) )-->
aRNAque is a simple, but efficient evolutionary algorithm for inverse RNA folding inspired by Lévy flights
For a given target structure in a dot-bracket representation, the tool allows to generate good quality (low ED and MFE) RNA sequences with the corresponding strucure close to the input target. The method relies on local mutations of nucleotide and base pairs independently with respect to some probabilities: P_N and P_C. More details about the choice of P_N and P_C is provide the SI of our paper.

The repo is organised as follows:
- [data](data/): The clean data used to produce the different plots presented in our paper.  The cleaned data are obtained by cleaning up the data generated from out file. for more details please refer to the python notebook [here](notebook/clean_data.ipynb)
- [docs] (docs): The files describing the evolutionary algorithm implemented.
- [images](images/): The plots (in pdf) used in the paper and the Python notebook code is in: [notebook/plots.ipynb](notebook/plots.ipynb).
- [src](src/): The source codes organized in three main parts:

    - [Utilities](src/utilities/): set of basic python functions usefull for our EA implementation and the script implementing the folding tool wrappers.
    
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
  - IPknot: the version we used in this work is the 0.0.5 and it can be downloaded to in the official website. We can also share with you on request the copy we have. 
  - HotKnots: the patched version of hotknots we used for our benchmark result can be found in the [thirdparty folder](thirdparty/Hotknots_v2.0_patched.zip). 
After installing the folding tools, make sure you have set an environmental variable for each tool to the bin decrectories. HOTKNOTS_ROOT for hotknots should be set to <path to hotknots>/bin  and IPKNOT for IPknot should be set to <path to the ipknot>/build.

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
                        values are {GC,GC1,GC2,GC3,GC4,GC25, GC50, GC75,ALL},
                        please check the online doc for more details (default:
                        GC2)
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
                        ip for IPknot and hk for Hotknots (default: v)
  --log                 Store the population for each instance of the inverse
                        folding in a folder (default: False)
  --verbose             Print the mean fitness evolution on a standard output
                        (default: False)
  --turner1999          Use the old energy parameters (default: False)
  -seed SEED            Seed for the initial population (default: None)
```
      $ python aRNAque.py --help

For pseudoknotted target, please choose the appropriate folding tool using the option -ft "ip" or -ft "hk".
