SHELL       = /bin/sh

# default config:
#    release
#    x86
#    shared
#    
#    to build 64 bit:                make -e M64=-m64
#    to build static:                make -e STATIC_LINK=1
#    to build 64 bit static debug:   make -e M64=-m64 STATIC_LINK=1 VARIANT=debug
#

VARIANT      = release
ifeq "$(VARIANT)" "debug"
    GCC_DEBUG_FLAGS += -g -D_DEBUG
endif

M64 = -m32
ARCH = x86
ifeq "$(M64)" "-m64"
    ARCH=x64
endif
  
EXAMPLE     = Normals
BINDIR      = ../bin/$(ARCH)/$(VARIANT)/$(EXAMPLE)
OBJDIR      = ../obj/$(ARCH)/$(VARIANT)/$(EXAMPLE)
LIBDIR      = ../lib/gcc4/$(ARCH)/$(VARIANT)
INCDIR      = ../include
COMMONDIR   = ../Common

TARGET      = $(BINDIR)/$(EXAMPLE)
OBJS        = main.o $(COMMONDIR)/Common.o
OBJSFINAL   = $(OBJDIR)/main.o $(OBJDIR)/Common.o

CC          = gcc
LD          = gcc 
RM          = /bin/rm -f

CXXFLAGS    = $(M64) $(GCC_DEBUG_FLAGS)
LDFLAGS     = $(M64) $(GCC_LINK_FLAGS)

RPATH_DEF    =
STATIC_LINK  =
FBXSDK_LIB   = -lfbxsdk
ifeq "$(STATIC_LINK)" ""
    CXXFLAGS += -DFBXSDK_SHARED
    RPATH_DEF = -Wl,-rpath $(PWD)/$(LIBDIR)
else
    FBXSDK_LIB = $(LIBDIR)/libfbxsdk.a
endif

LIBS = $(FBXSDK_LIB) -lm -lrt -luuid -lstdc++ -lpthread -ldl

.SUFFIXES: .cxx

all: $(TARGET)

$(TARGET): $(OBJS)
	mkdir -p $(BINDIR)
	$(LD) $(LDFLAGS) -o $@ $(OBJSFINAL) -L$(LIBDIR) $(LIBS) $(RPATH_DEF)
	cp *.fbx $(BINDIR)

$(COMMONDIR)/Common.o: $(COMMONDIR)/Common.cxx
	mkdir -p $(OBJDIR)
	$(CC) $(CXXFLAGS) -I$(INCDIR) -c $(COMMONDIR)/Common.cxx -o $(OBJDIR)/Common.o

.cxx.o:
	mkdir -p $(OBJDIR)
	$(CC) $(CXXFLAGS) -I$(INCDIR) -c $< -o $*.o
	mv $*.o $(OBJDIR)

clean:
	rm -rf $(BINDIR)
	rm -rf $(OBJDIR)

