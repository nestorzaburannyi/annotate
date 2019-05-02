from myoperator import RowOperator
import sys

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical.
class genus(RowOperator) :
    """Extract taxon from lineage to taxon_depth. Depth set to -2 if Archaea or Bacteria, depth 10 if Mammalia, 7 if Viridiplantae, 4 if other Eukaryota"""
    def __init__(self,glob):
        # define nicknames for column indices
        [self.species_col,self.lineage_col, self.taxon_col,self.kingdom_col]=glob.use_sheet("data").use_columns(["species","lineage","taxon","kingdom"])

    def process(self, row):
        x=row[self.lineage_col]
        #print >> sys.stderr, '#taxon:',x,self.lineage_col,row
        tmp=x.split("; ")
        taxon_depth=3
        lentmp=len(tmp)
        if lentmp<1: return
        kingdom=tmp[0]
        row[self.kingdom_col]=kingdom
        if kingdom == "Eukaryota":
          try:
                if tmp[1] == "Viridiplantae":
                        taxon_depth=6
                elif tmp[6] == "Mammalia":
                        taxon_depth=9
                elif tmp[12] == "Aves":
                        taxon_depth=14
          except:
                pass
        else:
                taxon_depth=lentmp-1
        if lentmp > taxon_depth:
                row[self.taxon_col]=tmp[taxon_depth].replace("Candidatus ","")
        else:
                row[self.taxon_col]=tmp[-1]

