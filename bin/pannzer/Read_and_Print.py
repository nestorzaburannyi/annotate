from __future__ import print_function
#################################
# Collect data reading
# and printing functions here
#
# Code mainly for SansPanz-project
############################################

import sys,re


##########################
# Debug libraries below  #
# Remove from final code #
##########################
# from python_diagnostics import *

###########################################
# Collect commands for reading data here
# ....
###########################################


def read_dict_data( file_name, col1, col2, header = False,
                    UseColNames = False, split_char = "\t", DoStrip = True, DoUpper = False, verbose=False):
    """
    Code is developed for Dict type data (Key = Value pairs)
    Usage DictOut = read_dict_data(file_name, col1, col2,
    ...                            header = True/False, UseColNames = True/False)
    input parameters:
    file_name
    col1        Column for key values (number or column name {UseColNames = True})
    col2        Column for values     (number or column name)
    header      True/False Does the file have a header
    UseColNames True/False Is the column selected with its name
    ...         Note that header must be True for UseColNames to be useful
    DoStrip     True/False Remove extra spaces (start and end) from value and key strings
    DoUpper     True/False Convert dictionary key to uppercase

    Stripping is not implemented for header row processing
    Bug generator: Petri Toronen
    """
    #################
    # Sub-functions #
    #################
    def num(str_in):
      # Code from Stack overflow
      # http://stackoverflow.com/questions/379906/parse-string-to-float-or-int
      # used by read_dict_data
      try:
         return int(str_in)
      except ValueError:
         try:
           return float(str_in)
         except ValueError:
           print("# Col value %s not a number" % str_in)
           return
    def FindHeaderMatch(header_row, find_str):
       # Test if header_row has a find_str
       # Separate function for processing exceptions
       # return the index
       try:
         pos1 = header_row.index(find_str)
       except ValueError:
         print("# These caused error:")
         print('# ',header_row)
         print("# \"%s\"" % find_str)
         print("# x not in list")
         pos1 = None
       return pos1
    #######################
    # Sub-functions ended #
    #######################

    if verbose: print("# read_dict_data",file_name,col1,col2,header,UseColNames,split_char,DoStrip,DoUpper)

    try:
      FILE_READ = open(file_name, "r")
    except IOError:
      print('# '+file_name +" was not found")
      print("# ReadDictData")
      return None       #THis will be the error case

    Data_dict = {};
    counter = 1      # This is used to filter header rows
    header_row = []
    if header:
      data_rows = 0
    else:
      data_rows = 1

    with open(file_name,"r") as FILE_READ:
     while 1:
       line=FILE_READ.readline()
       if not line: break
       if(data_rows == 0 & header):  # Process the header (if needed)
         header_row = line.split(split_char)
         header_row = list(map(str.strip, header_row)) #OR header_row = [x.strip(' ') for x in header_row]
       else:
         if DoStrip:
                tmp_table=line.strip().split(split_char)
         else:
                tmp_table  = line.split(split_char)
         if verbose and data_rows<10: print('# ',tmp_table)
         if(header and UseColNames and data_rows == 1):
           pos1 = FindHeaderMatch(header_row, col1)
           pos2 = FindHeaderMatch(header_row, col2)
           if verbose: print("# pos1=",pos1," pos2=",pos2)
         elif(data_rows == 1): # counter == 1
           pos1 = num(col1)
           pos2 = num(col2)
           if verbose: print("# pos1=",pos1," pos2=", pos2)

         #IF ELSE HEADER ends
         aa = tmp_table[pos1]
         try:
                bb = tmp_table[pos2]
         except:
                bb = aa # hack: map to self
         if DoUpper: aa=aa.upper()
         try:
          if DoStrip:
           Data_dict[aa.strip()] = bb.strip()
          else:
           Data_dict[aa] = bb
         except:
          print("# ERROR",tmp_table)
         if verbose and counter<10: print("# key = %s value = %s" %(aa,bb))
       # IF ELSE
       counter += 1
       data_rows += 1
    # FOR
    return Data_dict

def read_dict_counts(file_name, col1, col2):
    # special function to load word counts from text file (no tabs)
    counts={}
    print("# read_dict_count: %s" %file_name)
    with open(file_name,"r") as fh:
         while 1:
                line=fh.readline().lstrip()
                if not line: break
                data=line.split(" ",1)
                counts[data[col1].strip().upper()]=data[col2]
    return counts

def read_dict_GOdict(file_name):
    # special function to read GO dictionary: 3 tab-separarated columns, value is columns 2+3 (annotations + parents)
    GOdict={}
    with open(file_name,"r") as fh:
        while 1:
                line=fh.readline().strip()
                if not line: break
                row=line.split("\t",1)
                accession=row[0].strip()
                GOdict[accession]=row[1]
    return GOdict

def read_dict_PHR(file_name):
    # special function to read uniprot.phr: key is parsed from db|acc|id
    PHR={}
    with open(file_name,"r") as fh:
        while 1:
                line=fh.readline().strip()
                if not line: break
                row=line.split(" ",1)

                try:
                        (db,accession,pid)=row[0].split('|')
                        #LH remove uniprot tags OS=... GN=... PE=...
                        description=re.sub(r' \w{2}=.*$',"",row[1])
                        PHR[accession]=description
                        PHR[pid]=description
                except:
                        sys.stderr.write("PHR error %s\n" %line)
    return PHR


###########################################
# Folder-path-filename cleaning function
###########################################
def Remove_extra_slash(string_in):
   string_out = re.sub(r'/+', '/', string_in)
   return string_out

########################################################
# Now comes various random old data reading functions  #
########################################################


def read_ID_and_GO(file_name, name_col, GO_col,
                   header_rows = 0, split_char = "\t"):
  """
  This opens the file and read GO classes and Gene Identifier
  Everything is stored as a dict that has Gene ID as key
  GO classes are as a list (or array, how is it called...)
  Bug Generator: Petri Toronen

  read_ID_and_GO(file_name, name_col, GO_col,
  ...             header_rows = 0, split_char = "\t")
  """

  FILE_READ = open(file_name, "r")
  ID_GO_dict = {};
  counter = 0      # This is used to filter header rows

  for line in FILE_READ:
     if(counter >= header_rows):
       split_array = line.split(split_char)
       #print split_array
       if len(split_array) > GO_col:
          GO_annos = re.split(r'[ ,;.]+', split_array[GO_col].strip())
       else:
          GO_annos = [];
       for l in range(len(GO_annos)):
           GO_annos[l] = re.sub("GO:", "", GO_annos[l])
           GO_annos[l] = re.sub(r'\s+', "", GO_annos[l])
       seq_ID = re.sub(r'\s+','', split_array[name_col])
       ID_GO_dict[seq_ID] = GO_annos
       #print ID_GO_dict[seq_ID]
     #IF COUNTER
     counter += 1
  # FOR LINE IN .. ENDED
  # print_variables2(vars())
  FILE_READ.close()
  return  ID_GO_dict
# read_seqs_and_GO ended

def read_ID_and_GO_multiCols(file_name, name_col, GO_col,
                             header_rows = 0, split_char = "\t"):
  """
  This opens the file and read GO classes and Gene Identifier
  Everything is stored as a dict that has Gene ID as key
  GO classes are as a list (or array, how is it called...)
  Bug Generator: Petri Toronen
  THIS MODIFIED FROM THE OTHER CODE. HERE GO_col is a list of length 2
  """

  try:
    FILE_READ = open(file_name, "r")
  except IOError:
    print(file_name + " was not found")
    return None
  ID_GO_dict = {};
  counter = 0      # This is used to filter header rows

  for line in FILE_READ:
     if(counter >= header_rows):
       split_array = line.split(split_char)
       #print split_array
       if len(split_array) > GO_col[0]:    #FIRST COLUMN
          GO_annos1 = re.split(r'[ ,;.]+', split_array[GO_col[0]].strip())
       else:
          GO_annos1 = [];
       if len(split_array) > GO_col[1]:    #SECOND COLUMN
          GO_annos2 = re.split(r'[ ,;.]+', split_array[GO_col[1]].strip())
       else:
          GO_annos2 = [];
       GO_annos = list(set( GO_annos1 + GO_annos2 ))  # Join lists. REmove repetitive IDs
       GO_annos = list(set(GO_annos))
       for l in range(len(GO_annos)):
           GO_annos[l] = re.sub("GO:", "", GO_annos[l])
           GO_annos[l] = re.sub(r'\s+', "", GO_annos[l])
       seq_ID = re.sub(r'\s+','', split_array[name_col])
       ID_GO_dict[seq_ID] = GO_annos
       #print ID_GO_dict[seq_ID]
     #IF COUNTER
     counter =+ 1
  # FOR LINE IN .. ENDED
  # print_variables2(vars())
  FILE_READ.close()
  return  ID_GO_dict
# read_seqs_and_GO ended

def read_GO_ID_and_Par(file_name, name_col, parent_col,
                       header_rows = 0, split_char = "\t"):

  """
  This opens the file and read GO class and list of its parents
  Everything is stored as a dict that has GO class ID as a list

  """

  FILE_READ = open(file_name, "r")
  ID_GO_dict = {};
  counter = 0

  for line in FILE_READ:
     if counter >= header_rows:
       split_array = line.split(split_char)
       #print split_array
       if len(split_array) > parent_col:
          GO_par_list = re.split(r'[ ,;.]+', split_array[parent_col].strip())
       else:
          GO_annos = [];
       for l in range(len(GO_par_list)):
          GO_par_list[l] = re.sub("GO:", "", GO_par_list[l])
          GO_par_list[l] = re.sub(r'\s+', "", GO_par_list[l])
       GO_ID = re.sub("GO:",'', split_array[name_col])
       GO_ID = re.sub(r'\s+','', GO_ID)
       GO_par_list.append( GO_ID )
       ID_GO_dict[GO_ID] = GO_par_list
       #print ID_GO_dict[seq_ID]
     #IF COUNTER ... ENDED
     counter =+ 1
  # FOR LINE IN .. ENDED
  # print_variables2(vars())
  FILE_READ.close()
  return  ID_GO_dict
# read_seqs_and_GO ended


def read_ID_and_GO2(file_name, name_col, GO_col,
                    header_rows = 1, split_char = "\t"):

  """
  This opens the file and read GO classes and Gene Identifier
  This stores everything as array
  """

  FILE_READ = open(file_name, "r")
  ID_GO_list = [];
  counter = 0

  for line in FILE_READ:
     if counter >= header_rows:
       split_array = line.split(split_char)
       #print split_array
       if len(split_array) > GO_col:
          GO_annos = re.split(r'[ ,;.]+', split_array[GO_col])
       else:
          GO_annos = [];
       GO_annos[0] = re.sub("GO:", "", GO_annos[0])
       GO_annos[0] = re.sub(r'\s+', "", GO_annos[0])
       seq_ID = re.sub(r'\s+','', split_array[name_col])
       ID_GO_list.append( [ seq_ID, GO_annos[0] ] )
     # IF COUNTER < HEADER_ROWS
     counter += 1
     #print ID_GO_dict[seq_ID]
  # FOR LINE IN .. ENDED
  # print_variables2(vars())
  FILE_READ.close()
  return  ID_GO_list
# read_seqs_and_GO ended
