APP=hfr

FC=nvfortran
FC_FLAGS=-module ${MODDIR} -cuda -Minfo -cpp ${FC_OPT_FLAGS}

MODDIR=mod
SRCDIR=src
OUTDIR=out

all: samplemod
	${FC} ${FC_FLAGS} ${SRCDIR}/main.f08 -o ${OUTDIR}/${APP}

run:
	./${OUTDIR}/${APP}

samplemod: ${SRCDIR}/samplemod.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/samplemod.f08 -o ${OUTDIR}/samplemod.o

pgto: ${SRCDIR}/pgto.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/pgto.f08 -o ${OUTDIR}/pgto.o

sto_ng: pgto ${SRCDIR}/sto_ng.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/sto_ng.f08 -o ${OUTDIR}/sto_ng.o
