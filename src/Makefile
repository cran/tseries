LIBNAME = tsutils
LD = ld
CFLAGS += -Wall

OBJS = bdstest.o muin2ser.o tsutils.o

$(LIBNAME).so: $(OBJS)
	$(LD) $(SHLIBLDFLAGS) -o $@ $(OBJS) -lm 

clean:
	@rm -f *.o *.so

realclean: clean 
	@rm -f *~
