##
RCKSKEL_ROOT := ../rckskel
RCCEROOT = $(RCKSKEL_ROOT)/RCCE_V1.0.13
include $(RCCEROOT)/common/symbols
RCCE_ARCHIVE := $(RCCEROOT)/bin/SCC_LINUX/libRCCE_bigflags_nongory_nopwrmgmt.a
RCKSKEL_ARCHIVE := $(RCKSKEL_ROOT)/src/rckskel.a
RCKSKEL_INC := $(RCKSKEL_ROOT)/src
##
SOURCES := 	mcpsc.C \
			ce_align.C \
			pom/pdbutil.C  \
			pom/pom.C \
			pom/ipdb.C \
			cmp_util.C \
			pom/miscutil.C \
			pom/jdate.C \
			pom/ipdb_df.C \
			pom/linkedid.C \
			pom/derive.C \
			tmalign.C \
			usm.C

MODULES := 	ce_align.C \
			pom/pdbutil.C  \
			pom/pom.C \
			pom/ipdb.C \
			cmp_util.C \
			pom/miscutil.C \
			pom/jdate.C \
			pom/ipdb_df.C \
			pom/linkedid.C \
			pom/derive.C \
			tmalign.C \
			usm.C

PROGRAM := mcpsc
OBJECTS := $(SOURCES:.C=.o)
MOBJECTS := $(MODULES:.C=.o)

CC := g++
CE_CFLAGS := $(CFLAGS) -DFUNCPROTO -DREAD_WRITE -DSGI -O3 -std=c++11
CE_LFLAGS := $(LFLAGS) -lm -ffast-math -lboost_iostreams -lboost_system -lpthread -lz -static 

pc_test_all_v_all: $(MOBJECTS) 
	$(CC) $(CE_CFLAGS) -Wall -I$(RCKSKEL_INC) -c mcpsc.C -o mcpsc.o -D PC_TEST 
	$(CC) -o mcpsc mcpsc.o $(MOBJECTS) $(CE_LFLAGS) 

pc_test_usm: $(MOBJECTS) 
	$(CC) $(CE_CFLAGS) -Wall -c mcpsc.C -o mcpsc.o -D PC_TEST -D ONLY_USM 
	$(CC) -o mcpsc mcpsc.o $(MOBJECTS) $(CE_LFLAGS) 
	
pc_test_tmalign: $(MOBJECTS) 
	$(CC) $(CE_CFLAGS) -Wall -I$(RCKSKEL_INC) -c mcpsc.C -o mcpsc.o -D PC_TEST -D ONLY_TMALIGN
	$(CC) -o mcpsc mcpsc.o $(MOBJECTS) $(CE_LFLAGS) 
	
pc_test_ce: $(MOBJECTS) $(RCKSKEL_ARCHIVE)
	$(CC) $(CE_CFLAGS) -Wall -I$(RCKSKEL_INC) -c mcpsc.C -o mcpsc.o -D PC_TEST -D ONLY_CE
	$(CC) -o mcpsc mcpsc.o $(MOBJECTS) $(CE_LFLAGS) 	

pc_test_many_v_many: $(MOBJECTS) 
	$(CC) $(CE_CFLAGS) -Wall -I$(RCKSKEL_INC) -c mcpsc.C -o mcpsc.o -D PC_TEST -D MANY_MANY
	$(CC) -o mcpsc mcpsc.o $(MOBJECTS) $(CE_LFLAGS)
	
scc_1_test_many_v_many: $(MOBJECTS) $(RCKSKEL_ARCHIVE) $(RCCE_ARCHIVE)
	$(CC) $(CE_CFLAGS) -Wall -I$(RCKSKEL_INC) -c mcpsc.C -o mcpsc.o -D PC_TEST -D MANY_MANY -D SCC_1_CORE_TEST
	$(CC) -o mcpsc mcpsc.o $(MOBJECTS) $(RCKSKEL_ARCHIVE) $(RCCE_ARCHIVE) $(CE_LFLAGS)
	                
scc_test_all_v_all: $(MOBJECTS) $(RCKSKEL_ARCHIVE) $(RCCE_ARCHIVE)
	$(CC) $(CE_CFLAGS) -Wall -I$(RCKSKEL_INC) -c mcpsc.C -o mcpsc.o -D GREEDY_SPLIT
	$(CC) -o mcpsc mcpsc.o $(MOBJECTS) $(RCKSKEL_ARCHIVE) $(RCCE_ARCHIVE) $(CE_LFLAGS)

scc_test_many_v_many: $(MOBJECTS) $(RCKSKEL_ARCHIVE) $(RCCE_ARCHIVE)
	$(CC) $(CE_CFLAGS) -Wall -I$(RCKSKEL_INC) -c mcpsc.C -o mcpsc.o -D MANY_MANY -D GREEDY_SPLIT
	$(CC) -o mcpsc mcpsc.o $(MOBJECTS) $(RCKSKEL_ARCHIVE) $(RCCE_ARCHIVE) $(CE_LFLAGS)
	                
$(PROGRAM): $(OBJECTS) $(RCKSKEL_ARCHIVE) $(RCCE_ARCHIVE)
	@echo [$(OBJECTS)] are the objects
	$(CC) -o $@ $(OBJECTS) $(RCKSKEL_ARCHIVE) $(RCCE_ARCHIVE) $(CE_LFLAGS)

clean:
	rm -f *.o pom/*.o core $(PROGRAM)

.C.o: $(RCCEINCLUDE)/RCCE.h $(RCKSKEL_INC)/rckskel.h
	$(CC) $(CE_CFLAGS) -Wall -I$(RCKSKEL_INC) -c $< -o $@

##
