###################################################
#
# Makefile for pds_dmm
# Creator [Xcode -> Makefile Ver: Feb 14 2007 09:18:41]
# Created: [Mon Apr 14 08:04:30 2014]
#
###################################################

#
# Macros
#

CC = /usr/bin/g++
CC_OPTIONS = -O3
LNK_OPTIONS = 


#
# INCLUDE directories for pds_dmm
#

INCLUDE = -I.\
		-Ipds_dmm


#
# Build pds_dmm
#

pds_dmm : \
		./pds_dmm.o\
		./qFinderDMM.o\
		./linearalgebra.o
	$(CC) $(LNK_OPTIONS) \
		./pds_dmm.o\
		./qFinderDMM.o\
		./linearalgebra.o\
		-o pds_dmm

clean : 
		rm \
		./pds_dmm.o\
		./qFinderDMM.o\
		./linearalgebra.o\
		pds_dmm

install : pds_dmm
		cp pds_dmm pds_dmm

#
# Build the parts of pds_dmm
#


# Item # 1 -- pds_dmm --
./pds_dmm.o : pds_dmm.cpp
	$(CC) $(CC_OPTIONS) pds_dmm.cpp -c $(INCLUDE) -o ./pds_dmm.o


# Item # 2 -- qFinderDMM --
./qFinderDMM.o : qFinderDMM.cpp
	$(CC) $(CC_OPTIONS) qFinderDMM.cpp -c $(INCLUDE) -o ./qFinderDMM.o


# Item # 3 -- linearalgebra --
./linearalgebra.o : linearalgebra.cpp
	$(CC) $(CC_OPTIONS) linearalgebra.cpp -c $(INCLUDE) -o ./linearalgebra.o


##### END RUN ####
