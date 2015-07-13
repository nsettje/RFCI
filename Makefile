ICC := /opt/intel/composer_xe_2011_sp1.6.233/bin/intel64/icc
ICCLIB := /opt/intel/composer_xe_2011_sp1.6.233/compiler/lib/intel64
MKLINC := /opt/intel/composer_xe_2011_sp1.6.233/mkl/include
MKLROOT := /opt/intel/composer_xe_2011_sp1.6.233/mkl/
MKLLIB := /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64
INCLUDES := -isystem $(MKLINC) -I/home/settje/RFCI/src/

LDFLAGS :=-Wl,--start-group -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -Wl,--end-group -lpthread -lm
	#-Wl,--start-group $(MKLLIB)/libmkl_intel_lp64.a\
        $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a\
        -Wl,--end-group $(ICCLIB)/libiomp5.a -lpthread

CFLAGS := -o 

####################################################################
#.PHONY : default
default : bins method lib

###	Build util archive
#UTILSRC=$(shell find src/util -name "*.cc" | sed 's/^src[/]//')
#UTILOBJ=$(UTILSRC:%.cc=obj/%.o)
#UTILDEP=$(UTILOBJ:%.o=%.d)
#-include $(UTILDEP)
#.PHONY : util
#util : lib/util.a
#lib/util.a: $(UTILOBJ)
#	@mkdir -p lib
#	@rm -f $@
#	@echo Building archive $@ from $^
#	@ar rcs $@ $^

###	Build lib archive
LIBSRC=$(shell find src/lib -name "*.cc" | sed 's/^src[/]//')
LIBOBJ=$(LIBSRC:%.cc=obj/%.o)
LIBDEP=$(LIBOBJ:obj/%.o=%.d)
-include $(LIBDEP)
.PHONY : lib
lib : lib/lib.a
lib/lib.a: $(LIBOBJ)
	@mkdir -p lib
	@rm -f $@
	@echo Building archive $@ from $^
	@ar rcs $@ $^

###	Build method archive
METHSRC=$(shell find src/method -name "*.cc" | sed 's/^src[/]//')
METHOBJ=$(METHSRC:%.cc=obj/%.o)
METHDEP=$(METHOBJ:obj/%.o=%.d)
-include $(METHDEP)
.PHONY : method 
method : lib/method.a
lib/method.a: $(METHOBJ)
	@mkdir -p lib
	@rm -f $@
	@echo Building archive $@ from $^
	@ar rcs $@ $^


.PHONY: bins
bins: RFCI
	@echo Building executable RFCI
	@rm -f *.d

RFCI: obj/core/CI.o lib/lib.a 
	@echo Linking $^
	@$(ICC) $(CFLAGS) $@ $(INCLUDES) $^ $(LDFLAGS) -openmp -I${MKLROOT}/include

###	Generic make rules
obj/%.o: src/%.cc
	@echo 1Compiling $^ into $@
	@mkdir -p $(dir $@)
	@$(ICC)     -c $(INCLUDES) $< -o $@
	@$(ICC) -M -c $(INCLUDES) $< > obj/$*.d.tmp
	@sed "0,/^.*:/s//$(subst /,\/,$@):/" obj/$*.d.tmp > obj/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < obj/$*.d.tmp | fmt -1 | \
		sed -e 's/^ *//' -e 's/$$/:/' >> obj/$*.d
	@rm -f obj/$*.d.tmp

###	Generic make rules
obj/%.o: src/%/*.cc 
	@echo 2Compiling $^ into $@
	@mkdir -p $(dir $@)
	@$(ICC)     -c $(INCLUDES) $< -o $@
	@$(ICC) -M -c $(INCLUDES) $< > $*.d.tmp
	@sed "0,/^.*:/s//$(subst /,\/,$@):/" $*.d.tmp > obj/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
		sed -e 's/^ *//' -e 's/$$/:/' >> obj/$*.d
	@rm -f $*.d.tmp


.PHONY : clean
clean :
	rm -f src/core/*.d src/core/*.o src/util/*.a
	rm -rf obj lib

