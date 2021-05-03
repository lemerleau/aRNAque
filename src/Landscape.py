

"""

"""
import RNA





class Landscape (object) :

    def __init__ (self, target) :
        self.target = target

    def fitness (self, structure) :
        if '{' in structure :
            str_ = structure.replace('{','(')
            str_ = str_.replace('}',')')
            return 1./(1+RNA.hamming_distance(self.target, str_))

        return 1./(1+RNA.hamming_distance(self.target, structure))

    def ens_defect(self, sequence) :

        fc = RNA.fold_compound(sequence)
        fc.pf()
        fc.bpp()
        ed = fc.ensemble_defect(self.target)
        return ed
