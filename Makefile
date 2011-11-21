include Makefile.in.mac

.PHONY:	all
#all:	matmul matmul-basic matmul-blocked matmul-blas matmul-pblas matmul-f2c
all:	matmul matmul-basic matmul-blocked matmul-blas matmul-pblas 

# ---
# Rules to build the drivers

matmul: $(OBJS) dgemm.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

matmul-basic: $(OBJS) basic_dgemm.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

matmul-blocked:	$(OBJS) blocked_dgemm.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

matmul-blas:	$(OBJS) blas_dgemm.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS) $(LIBBLAS)

matmul-pblas:	$(OBJS) blas_dgemm.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS) $(LIBPBLAS)

matmul-f2c:	$(OBJS) f2c_dgemm.o fdgemm.o
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS) 

# --
# Rules to build object files

%.o:%.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $<

%.o:%.f
	$(FC) -c $(FFLAGS) $<

blas_dgemm.o: blas_dgemm.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCBLAS) $< 

# ---
# This is just a suggestion on how to generate timing plots...  Feel
# free to improve on these, so long as you show MFlop/s v. matrix size.
# Submit matmul via qsub to generate the raw timing data.

.PHONY: run
run: matmul
	./make_sge.sh ./matmul mine
	qsub serial.qsub

.PHONY: run-basic
run-basic: matmul-basic
	./make_sge.sh ./matmul-basic basic
	qsub serial.qsub

.PHONY: run-blocked
run-blocked: matmul-blocked
	./make_sge.sh ./matmul-blocked blocked
	qsub serial.qsub

.PHONY: run-blas
run-blas: matmul-blas
	./make_sge.sh ./matmul-blas atlas
	qsub serial.qsub

.PHONY: run-pblas
run-pblas: matmul-pblas
	./make_sge.sh ./matmul-pblas patlas
	qsub serial.qsub

.PHONY: run-f2c
run-f2c: matmul-f2c
	./make_sge.sh ./matmul-f2c f2c
	qsub serial.qsub

# --
# Rules to compile timepgf (if everything is available)
timepgf.pdf: timepgf.tex
	pdflatex timepgf.tex
	pdflatex timepgf.tex

# ---

.PHONY:	clean realclean tgz
clean:
	rm -f matmul matmul-basic matmul-blocked matmul-blas matmul-pblas \
		matmul-f2c *.o
	rm -f serial.qsub*
	rm -f timepgf.aux timepgf.log timepgf.out

realclean:	clean
	rm -f *~ timing-*.dat timing-*.out timing-*.pdf

tgz: realclean
	(cd ..; tar -czf matmul.tgz matmul)
