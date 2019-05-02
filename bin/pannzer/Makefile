# Makefile for  SANS-server.
# # Copyright (C) 2014 Liisa Holm.
# 
# ##############################################################################################
# # start of configuration section                                                             #
# ##############################################################################################
#
 SHELL = /bin/sh
# Directories <--- modify these!!
 HOME=./
 SRC=$(HOME)
 BIN=$(HOME)
 DATA=/data/uniprot/
#
# # Data for experiments
 UNIPROT_TREMBL=${DATA}/uniprot_trembl.fasta.gz # uniprot_trembl sequences in FASTA format
 UNIPROT_SPROT=${DATA}/uniprot_sprot.fasta.gz # uniprot_swissprot sequences in FASTA format

# ##############################################################################################
# # end of configuration section                                                               #
# ##############################################################################################
#
# # Compiler options
 CC = gcc
 FC = gfortran
 CC_MPI = mpicc                   # C compiler
 FC_MPI = mpif90              # Fortran-95 compiler
 OPT = -O3                  # uncomment this for optimization
 DEB = # -g                 # uncomment this for debugging
 DEFS =                     # place other -D types of things here
 CFLAGS = $(OPT) $(DEB) $(DEFS)
 LFLAGS =
#
# # objects
 OBJS = server_c.o sais.o sais8.o 
# modules
 MODS =  sans.mod myhash.mod keyspace.mod maxheap.mod sanslocaldata.mod altschul.mod
#
 EXES = sais8 isa8 sap8 client server
#
#
# # Top-level targets
all: $(OBJS) $(EXES) 
install: all

#-----------------------------
# sans program suite
# #-----------------------------
#
# # parameters module
sais.o: $(SRC)/sais.c $(SRC)/sais.h
	$(CC) $(CFLAGS) -c $(SRC)/sais.c

server_c.o: $(SRC)/server_c.c
	$(CC_MPI) $(CFLAGS) -c $(SRC)/server_c.c

# suffix array generation 
sais8: $(SRC)/sais8.f $(SRC)/sais8.c $(SRC)/sais8.h
	$(CC) $(CFLAGS) -c $(SRC)/sais8.c 
	$(FC) $(CFLAGS) $(SRC)/sais8.f sais8.o -o sais8

isa8: $(SRC)/isa8.f 
	$(FC) $(CFLAGS) $(SRC)/isa8.f -o isa8

sap8: $(SRC)/sap8.f 
	$(FC) $(CFLAGS) $(SRC)/sap8.f -o sap8

# client
client: $(SRC)/client.c
	$(CC) $(CFLAGS) $(SRC)/client.c -o client

# server
server: $(SRC)/server.f server_c.o sais.o 
	$(FC_MPI) $(SRC)/server.f sais.o server_c.o -o server

# clean
clean: 
	rm -rf $(OBJS) $(MODS) $(EXES)

 
