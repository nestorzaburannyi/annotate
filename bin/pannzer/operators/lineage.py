import sys
from myoperator import RowOperator

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical.
class lineage(RowOperator) :
    """Do a dictionary lookup for the taxonomic lineage of a species."""
    def __init__(self,glob):
        self.glob=glob
        # define nicknames for column indices
        [self.species_col, self.lineage_col]=self.glob.use_sheet("data").use_columns(["species","lineage"])
        # use online dictionary. Object handles in glob are hardcoded
        self.glob.use_online_dictionaries(["LINEAGE"])

    def process(self,row):
        species=row[self.species_col]
        row[self.lineage_col]=self.get_lineage(species)

    def get_lineage(self,species):
        ucspecies=species.upper() # lineage keys are UPPERCASE species
        if ucspecies in self.glob.lineage:
                tmp=self.glob.lineage[ucspecies]
                if tmp: return(tmp)
        if species != "unknown": sys.stderr.write("# Warning: unknown species ||%s||\n" %species)
        return("unclassified")

