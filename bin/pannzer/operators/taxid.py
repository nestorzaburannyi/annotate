from __future__ import print_function
from myoperator import RowOperator
from Read_and_Print import read_dict_data
import sys

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical.
class taxid(RowOperator) :
    """Extract taxon from lineage to taxon_depth."""
    def __init__(self,glob):
        self.glob=glob
        # define nicknames for column indices
        [self.species_col, self.taxid_col]=glob.use_sheet("data").use_columns(["species","taxid"])
        # load taxid from dictionary
        #fn=glob.param['DATA_DIR']+'/'+glob.param['DATA_TAXONOMY']
        #self.taxid=read_dict_data(fn,3,0) # Scientific name, Taxon
        glob.use_online_dictionaries(['TAXID'])

    def process(self, row):
        species=row[self.species_col].upper() # keys are uppercase
        try:
                row[self.taxid_col]=self.glob.taxid[species]
        except:
                print('# species not found in taxonomy:',species, file=sys.stderr)

