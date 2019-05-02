from myoperator import BlockOperator
import sys
from PannzerFunctions import createCorrectedTaxDistances

class Taxonomy(BlockOperator):
        """
        Calculate taxdist function between query and sbjct species. Assume SANS tabular output.

        Creates data column 'taxdist'.
        Inputs: data column 'species'

        Species comparisons in upper case
        """
        def __init__(self, glob):
                sys.stderr.write("# Init Taxonomy\n")
                self.glob=glob
                [self.taxdist_col, self.species_col, self.lineage_col]=glob.use_sheet("data").use_columns(["taxdist","species","lineage"])
                [self.a]=glob.use_operators(['lineage'])
                # locals
                self.taxdist={} # for hashing td of recurrent lineages relative to queryspecies
                self.taxdistfunction=createCorrectedTaxDistances()
                self.previous_queryspecies=''
                self.querytd=None

        def process(self,block):
                # add lineage col to hits
                for row in block: self.a.process(row) # add lineage column
                # use first unclassified lineage from hits list
                self.querylineage=""
                self.queryspecies=self.glob.QUERYSPECIES
                for row in block:
                        x=row[self.lineage_col]
                        if x =="unclassified": continue
                        self.querylineage=x
                        self.queryspecies=row[self.species_col]
                        #sys.stderr.write("# Query's best-guess-lineage is %s\n" %(self.querylineage))
                        break
                self.querytd=len(self.querylineage.split(";")) # length of vector
                if self.previous_queryspecies != self.queryspecies:
                        self.taxdist={}
                        self.previous_queryspecies=self.queryspecies
                for row in block: self.process_row(row)

        def process_row(self,row):
                line=row[self.lineage_col]
                if line in self.taxdist:
                        td=self.taxdist[line] # use cached value
                else:
                        pl=self.get_pathlength(line)
                        td=self.taxdistfunction[pl]
                        self.taxdist[line]=td
                row[self.taxdist_col]=str(td) # data is string

        def get_pathlength(self,lineage):
                if len(self.querylineage)<1: return(0)
                x=lineage.split(";")
                std=len(x)
                td=self.querytd+std
                i=0
                while i<self.querytd and i<std and self.querylineage[i] == x[i]:
                        td=td-2
                        i=i+1
                return(td)

