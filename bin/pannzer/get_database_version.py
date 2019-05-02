from __future__ import absolute_import
from __future__ import print_function
import Runner, Parameters,XMLParser
glob=Parameters.WorkSpace()
glob.param["input_OUTFILES"]=""
z=Runner.Runner(glob,operator_name="SANS")
z.xml=XMLParser.XMLParser()
metadata=z.SANSquery(">dummy\nHAHAHAHAHAHAHAHAHAHAHAHAHAHAHA\n")[1]
print("%s consisting of %s letters and %s sequences" %(metadata[0],metadata[1],metadata[2]))
