CC=g++
ROOTGEN=rootcint
CFLAGS=-g -Wall `root-config --cflags`
CPPFLAGS=-I ./include
LDFLAGS=`root-config --glibs`

ROOTINCLDIR=./
INCLDIR=./include
SRCDIR=./src
OBJDIR=./objs
BINDIR=./bin

SRC=$(wildcard $(SRCDIR)/*.cpp)
OBJS=$(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

DICTOBJ=$(OBJDIR)/kinematics_dict.o
DICTSRC=$(SRCDIR)/kinematics_dict.cxx
DICT_PAGES=$(INCLDIR)/Kinematics.h $(INCLDIR)/LinkDef_Kinematics.h

EXE=$(BINDIR)/mask 

CLEANUP=$(EXE) $(OBJS) $(DICTOBJ) $(DICTSRC)

.PHONY: all clean

all: $(EXE) $(TPEXE)

$(EXE): $(DICTOBJ) $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(DICTOBJ): $(DICTSRC)
	$(CC) $(CFLAGS) $(CPPFLAGS) -I $(ROOTINCLDIR) -c $^ -o $@
	mv $(SRCDIR)/*.pcm $(BINDIR)

$(DICTSRC): $(DICT_PAGES)
	$(ROOTGEN) -f $@ $^

VPATH= $(SRCDIR):./testplots/
$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	$(RM) $(CLEANUP) $(BINDIR)/%.pcm
