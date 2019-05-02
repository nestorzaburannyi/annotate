# TextOperator can be used to parse text files and output tab

from myoperator import TextOperator
import sys,re

class obo(TextOperator):
    """
Input: go-basic. obo
Output: goid, directParentSet, CoParentSets
    """
    def __init__(self,glob):
        # set parameter
        self.glob=glob
        [self.data,self.summary_data]=glob.use_sheets(['data','goidsets'])
        self.cols1=self.data.use_columns(['goid','mapped_to','children','directParentSet','CoParentSets','ontology','desc'])
        self.cols2=self.summary_data.use_columns(['goidset'])
        # do sums over input stream
        self.goid=None
        self.goids=[] # list of goids
        self.mapped_to={}
        self.children={}
        self.parents={} # key=goid
        self.coparents={} # key=goid
        self.p=[] # list of direct parents
        self.obsolete=True
        self.altids=[] # key=goid; altids have no obo-entry of their own
        self.namespace={}
        self.desc={}

    def process(self, line):
        if line[0:7]=='id: GO:':
                self.goid=line[7:14]
                self.p=[]
                self.obsolete=False
                self.altids=[]
                self.godesc=""
        elif line[0:9]=='is_a: GO:':
                self.p.append(line[9:16])
        elif line[0:25]=='relationship: part_of GO:':
                self.p.append(line[25:32])
        elif line[0:17]=='is_obsolete: true':
                self.obsolete=True
        elif line[0:11]=='alt_id: GO:':
                self.altids.append(line[11:18])
        elif line[0:6]=='name: ':
                self.godesc=line[6:].rstrip()
        elif line[0:11]=='namespace: ':
                x=line[11:].rstrip()
                if x=='biological_process':
                        self.ontology='BP'
                elif x=='molecular_function':
                        self.ontology='MF'
                elif x=='cellular_component':
                        self.ontology='CC'
                else:
                        self.ontology=''
        elif line=='\n':
                # term is complete
                if self.obsolete: return
                # copy data to alt_ids
                for x in [self.goid]+self.altids: self.savegoid(x,self.goid,self.p,[],self.ontology,self.godesc)

    def savegoid(self,goid,mapped_to,p,cop,ontology,godesc):
        if goid in self.goids: return # no overwriting!
        self.goids.append(goid)
        self.mapped_to[goid]=mapped_to
        p.sort()
        self.parents[goid]=p # list of parents where goid is child
        self.coparents[goid]=cop # list of coparents where goid is parent
        self.namespace[goid]=ontology
        self.desc[goid]=godesc
        for parentgoid in p:
                if not parentgoid in self.children: self.children[parentgoid]=[]
                self.children[parentgoid].append(goid)

    def getcoparentsets(self,parentset):
        for goid in parentset:
                cop={}
                for x in parentset:
                        cop[x]=1
                tmp=list(cop.keys())
                tmp.sort()
                self.coparents[goid].append(",".join(tmp))

    def finalise(self):
        for goid in self.goids:
                self.getcoparentsets(self.parents[goid])
        # collect unique goidsets for which we need counts
        u={}
        for goid in self.goids:
                u[goid]=1
                for x in self.coparents[goid]: u[x]=1 # stringified coparent set
        tmp=list(u.keys())
        tmp.sort()
        for x in tmp:
                self.summary_data.append_row([x])
        # remove redundancy from coparents
        for goid in self.goids:
                u={}
                for x in self.coparents[goid]: u[x]=1 # stringified coparent set
                tmp=list(u.keys())
                tmp.sort()
                self.coparents[goid]=tmp
        # write to speadsheet...
        for goid in self.goids:
                if not goid in self.children: self.children[goid]=[]
                self.data.append_row([goid,self.mapped_to[goid],",".join(self.children[goid]),",".join(self.parents[goid]),";".join(self.coparents[goid]),self.namespace[goid],self.desc[goid]])
        # stream summary output
        self.data.output(result=True)
        self.summary_data.output(result=True)

