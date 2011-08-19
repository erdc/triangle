CC         = gcc
LD         = gcc
SWIG       = swig
SRCDIR     = ../src
SRCS       = triangleModule.c
OBJS       = triangleModule.o
INCLUDES   = -I. -I${SRCDIR} 
TARGET     = triangulate
WRAP       = 
INTERFACE  = 
#if have multiple files with Array Api may need to define
#PY_ARRAY_UNIQUE_SYMBOL
DEFS       = -DDEBUG -DMWF_ADDED_FLAGS
PYTHON 	   = /usr/bin/python
PYINCLUDE  = /Library/Frameworks/Python.framework/Versions/2.4/include/python2.4
darwin:
	$(CC) -g ${INCLUDES} -I${PYINCLUDE} $(DEFS) -c  $(SRCS) 
	$(CC) -g ${INCLUDES} -I${PYINCLUDE} -bundle -bundle_loader $(PYTHON) -framework Python $(OBJS)  -lm -o $(TARGET).so -L../lib -ltri

wrappers:
	$(SWIG)  -I/usr/include -python ${INCLUDES} ${INTERFACE}

clean::
	rm -f $(WRAP).o $(TARGET).so $(TARGET).pyc *.o

clobber::clean
	rm -f $(TARGET).py  
