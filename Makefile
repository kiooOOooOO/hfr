APP=hfr

FC=nvfortran
FC_FLAGS=-stdpar=multicore -module ${MODDIR} -lblas -cpp ${FC_OPT_FLAGS}

MODDIR=mod
SRCDIR=src
OUTDIR=out

all: hf
	${FC} ${FC_FLAGS} ${SRCDIR}/main.f08 ${OUTDIR}/*.o -o ${OUTDIR}/${APP}

run:
	./${OUTDIR}/${APP}

clean:
	rm -f ${OUTDIR}/*.o ${OUTDIR}/${APP}

samplemod: ${SRCDIR}/samplemod.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/samplemod.f08 -o ${OUTDIR}/samplemod.o

pgto: ${SRCDIR}/pgto.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/pgto.f08 -o ${OUTDIR}/pgto.o

sto_ng: pgto ${SRCDIR}/sto_ng.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/sto_ng.f08 -o ${OUTDIR}/sto_ng.o

hf: sto_ng matrix ${SRCDIR}/hf.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/hf.f08 -o ${OUTDIR}/hf.o

matrix: ${SRCDIR}/matrix.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/matrix.f08 -o ${OUTDIR}/matrix.o
