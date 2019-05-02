from myoperator import BlockOperator
import sys

class BLAST2GO2(BlockOperator):
        """
        Implement blast2go
        """
        def __init__(self, glob,verbose=True):
                sys.stderr.write("# Init blast2go\n")
                self.glob=glob
                # two spreadsheets (input, output)
                [self.data,self.anno_data]=glob.use_sheets(["data","anno"])
                # define input column nicknames
                [self.qpid_col,self.status_col,self.bits_col,self.spid_col,self.golist_col,self.ecwlist_col,self.p_ident_col]= \
                     self.data.use_columns(['qpid','status','bits','spid','GOclass','evidence_code_weight','pide'])
                self.anno_data.use_columns(['qpid','type','score','PPV','id','desc'])
                if 'B2G_thresh' in glob.param:
                   self.threshold = float( glob.param['B2G_thresh'])  # allow float values also
                else: # use here default value
                   self.threshold = 55
                header_string = "Selected_classes_with_" + str(self.threshold)
                # sefine output column headers
#               self.goclass_data.use_columns(["qpid","goid","gocount","BLAST2GO_PPV","BLAST2GO_rank",
#                                              "BLAST2GO_score", "SelectedClasses", "BestBitScore"])
                self.MAXHITS=glob.param['PANZ_MAXHITS']
                # filter hit list; fetch go annotations and evidence-code weights
                [self.d,self.l]=glob.use_operators(['Filter','B2G'])
                # fetch annotations of go classes
                glob.use_online_dictionaries(["GOIDELIC"]) # to get ontology
                # method-specific parameter
                self.GO_weight=5
                #self.threshold = threshold


        def process(self,block):
                # sort SSRL
                if len(block)==0: return
                self.d.process(block)  # Filter status
                for row in block: self.l.process(row) # B2G is RowOperator; creates GOclass GOclass_count   evidence_code_weight columns

                # Added this loop
                n=0
                for row in block:   #PT This loop is copied from RM3 process
                        if row[self.status_col] == 'False': continue
                        row[self.status_col]="False"
                        if n>=self.MAXHITS: continue
                        if len(row[self.golist_col])>0:
                                row[self.status_col]="True"
                                n+=1
                                continue

                B2GO_dict1 = {}
                B2GO_dict2 = {}
                BestBlast_dict = {}      # This represents the best bit-score per GO class
                self.sizeInBackGround={} # GOclass_count per GOclass [could take directly from GOIDELIC]

                for row in block:
                   if row[self.status_col]=='False': continue
                   # split GOclass and evidence_code_weight lists to entries in goclass table
                   goidlist=row[self.golist_col].rstrip().split(' ')
                   goweights=row[self.ecwlist_col].rstrip().split(' ')
                   if len(goidlist) != len(goweights):
                        sys.stderr.write("# Error in BLAST2GO2.\n")
                        sys.stderr.write("# vectors goidlist and goweights don't match\n")
                        sys.stderr.write("# goidlist:\n# %s\n" %(", ".join(goidlist)))
                        sys.stderr.write("# goweights:\n# %s\n" %(", ".join(goweights)))
                        sys.stderr.write("# Current row: \n# %s\n" %("/ ".join(row)))
                   # Collect data required for each GO class
                   for goid in goidlist:
                      if not goid in self.sizeInBackGround:
                         try:
                            self.sizeInBackGround[goid]=int(self.glob.GOcounts[goid]) #int(gocntlist[i])
                         except KeyError:
                            sys.stderr.write("# no GOcounts for goid = %s\n" %(goid))
                            sys.stderr.write("# It will be removed from the analysis\n")
                            sys.stderr.write("# Current row: \n# %s\n" %("/ ".join(row)))
                         del_ind = goidlist.index(goid)
                         del goidlist[del_ind]  # Assume goid is obsolete. Better to remove it
                         del goweights[del_ind] # Must remove same weight item

                   seq_Ident = float(row[self.p_ident_col])
                   BitScore  = float(row[self.bits_col])
                   spid=row[self.spid_col]
                   qpid=row[self.qpid_col]
                   n=len(goidlist)
                   m=len(goweights)
                   if n > 0 and len(row[self.ecwlist_col]) > 0:   # PT: second check is a hack.

                       for i in range(n):

                                tmp_score = float(goweights[i])*seq_Ident*100
                                tmp_ID = goidlist[i]
                                if tmp_ID in B2GO_dict1:
                                   B2GO_dict2[tmp_ID] = B2GO_dict2[tmp_ID] + self.GO_weight

                                   if tmp_score > B2GO_dict1[tmp_ID]:
                                      B2GO_dict1[tmp_ID] = tmp_score
                                   if BitScore > BestBlast_dict[tmp_ID]:
                                      BestBlast_dict[tmp_ID] = BitScore
                                else:
                                   B2GO_dict1[tmp_ID] = tmp_score
                                   B2GO_dict2[tmp_ID] = 0 # Blast2GO uses (#GO - 1)
                                   BestBlast_dict[tmp_ID] = BitScore

                # Following steps select the
                # GO classes in blast2go style
                RemoveParents = {}      # These will be removed
                Accept_IDs = {}         # These are accepted results
                for go_id in B2GO_dict1:
                   if B2GO_dict1[go_id] + B2GO_dict2[go_id] > self.threshold \
                   and go_id not in RemoveParents:
                      Accept_IDs[go_id] = 1
                      if go_id in self.glob.GOparents:
                        parent_list = self.glob.GOparents[go_id]
                      else:
                        parent_list = []
                      for id2 in parent_list:
                        RemoveParents[id2] = 1
                        Accept_IDs.pop(id2, None)

                # write table with goid entries
                for go_id in B2GO_dict1:
                   # reject root classs
                   if self.sizeInBackGround[goid] == self.glob.rootcount[self.glob.ontology[goid]]: continue
                   #print B2GO_dict2[go_id], B2GO_dict1[go_id]
                   if not go_id in Accept_IDs: continue
                   datarow = [qpid, '%s_%s' %(self.glob.ontology[go_id],'B2GO'), str(B2GO_dict2[go_id]), '0.5', go_id,self.glob.godesc[go_id]]
                   #print datarow
                   self.anno_data.append_row(datarow)
