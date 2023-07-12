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
