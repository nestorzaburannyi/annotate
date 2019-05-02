from __future__ import print_function
#############################
# Data objects              #
#############################

import sys

class SpreadSheet:
        """
        A simple spreadsheet-like data object. Attributes:

                ncols - number of columns
                nrows - number of rows
                colnames - column headings. An array of ncols strings.
                block - an array of nrows data rows. A row is an array of ncols strings.
                fh - output channel. Default is sys.stdout
                row_status - True/False. Controls whether row is output.
                col_status - True/False. Controls whether column is output.

        All data elements are strings. Convert to appropriate type before arithmetics!
        """
        def __init__(self):
                # zero-by-zero matrix
                self.empty_header()
                self.empty_block()
                self.fh=sys.stdout # default output channel
                self.connection=None # used by socket-server
                self.pythonversion=list(sys.version_info)[0]

        def output(self,header=False,result=False):
                """
                Write column names to file handle or socket connection if header is True.
                Write block to file handle or socket connection if result is True.
                Hide column if col_status is False.
                Hide row if row_status is False.
                No output if file handle self.fh is undefined.

                """
                if self.connection or self.fh:
                  string=''
                  if self.pythonversion>2:
                        outcols=list(filter( lambda i: self.col_status[i], range(0,self.ncols) ))
                  else:
                        outcols=filter( lambda i: self.col_status[i], range(0,self.ncols) )
                  # build result string: header
                  if header:
                        string+=("\t".join( map( lambda a: self.colnames[a], outcols ) )+"\n")
                  # build result string: block
                  if result:
                        for i in range(0,self.nrows):
                                row=self.block[i]
                                # skip filtered row
                                if not self.row_status[i]: continue
                                # skip filtered columns
                                try:
                                  x = "\t".join( map( lambda a: row[a], outcols ) )+"\n"
                                  string += x
                                except:
                                  sys.stderr.write("# Warning from output: %s\n" %row)
                                  print(outcols)
                  # write to file handle
                  if self.connection:
                        self.connection.sendall(string)
                  else:
                        self.fh.write(string)

        def hide_from_output(self,hidecols):
                """
                Hide columns from output
                """
                for x in hidecols: self.col_status[self.colnames.index(x)]=False

        def empty_header(self):
                """
                Initialize colnames, col_status and ncols.
                """
                self.colnames=[] # new columns appended to input
                self.col_status=[] # Booleans, output column if True
                self.ncols=0

        def empty_block(self):
                """
                Initialize block, rownames, row_status and nrows.
                """
                self.block=[]
                self.rownames=[]
                self.row_status=[] # Booleans, output row if True
                self.nrows=0

        def sort_block(self,sortcol,reverse=False):
                """
                Numeric sort in place.
                """
                self.block.sort(key=lambda x: float(x[sortcol]), reverse=reverse)

        def append_row(self,row):
                """
                Create new row in block. Add data from row to the left and pad with empty elements.
                """
                self.block.append(row)
                self.row_status.append(True)
                for i in range(len(row),self.ncols): row.append('n.d.')
                self.nrows+=1

        def append_columns(self,newcolnames):
                """
                Grow existing block to the right with new columns and pad empty elements on each row.
                """
                n=self.create_columns(newcolnames)
                for row in self.block:
                        for i in range(0,n): row.append('n.d.')

        def create_columns(self,newcolnames):
                """
                Create new column for each element of newcolnames array. Save named index of column.
                Returns list of column indices of new columns created. Only modifies colnames (not block).
                """
                ix=[]
                for x in newcolnames:
                        self.colnames.append(x)
                        self.col_status.append(True)
                        ix.append(self.colnames.index(x))
                        self.ncols+=1
                return(ix)

        def use_columns(self,column_names):
                """Create columns if they don't exist yet. Return column indices."""
                ix=[]
                for x in column_names:
                        if not x in self.colnames:
                                self.colnames.append(x)
                                self.col_status.append(True)
                                self.ncols+=1
                        ix.append(self.colnames.index(x))
                return(ix)

        def get_col_index(self,column_names, verbose=False):
                """Returns list of column indices for list of column_names (None if column name not found)"""
                ix=[]
                for x in column_names:
                        try:
                                i=self.colnames.index(x)
                        except:
                                if verbose: sys.stderr.write("# Warning from get_col_index: column name %s not found\n" %x)
                                i=None
                        ix.append(i)
                return(ix)

