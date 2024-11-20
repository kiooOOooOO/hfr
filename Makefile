APP=hfr

FC=mpifort
FC_FLAGS=-module ${MODDIR} -lblas -cpp ${FC_OPT_FLAGS}
#FC_OPT_FLAGS?=-Dpure="" -g

MODDIR=mod
SRCDIR=src
OUTDIR=out

INFILE?=./inputs/h2.dat
DUMP?=no
OUTFILE?=result.dat

MPI_WORKER_NUM?=6

all: hf mulliken
	${FC} ${FC_FLAGS} ${SRCDIR}/main.f08 ${OUTDIR}/*.o -o ${OUTDIR}/${APP}

eritest: pgto
	${FC} ${FC_FLAGS} ${SRCDIR}/eritest.f08 ${OUTDIR}/*.o -o ${OUTDIR}/eritest

run:
	mpirun --mca btl_base_warn_component_unused 0 --mca mpi_cuda_support 0 -n ${MPI_WORKER_NUM} ./${OUTDIR}/${APP} ${INFILE} ${OUTFILE} ${DUMP} && \
	cat result.dat

run-mhost:
	mpirun --mca btl_base_warn_component_unused 0 --mca mpi_cuda_support 0 --hostfile hosts -n ${MPI_WORKER_NUM} ${OUTDIR}/${APP} ${INFILE} ${OUTFILE} ${DUMP} && \
		cat result.dat

clean:
	rm -f ${OUTDIR}/*.o ${OUTDIR}/${APP}

samplemod: ${SRCDIR}/samplemod.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/samplemod.f08 -o ${OUTDIR}/samplemod.o

pgto: ${SRCDIR}/pgto.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/pgto.f08 -o ${OUTDIR}/pgto.o

sto_ng: pgto ${SRCDIR}/sto_ng.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/sto_ng.f08 -o ${OUTDIR}/sto_ng.o

hf: hf_situation sto_ng matrix ${SRCDIR}/hf.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/hf.f08 -o ${OUTDIR}/hf.o

hf_situation: sto_ng timer ${SRCDIR}/hf_situation.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/hf_situation.f08 -o ${OUTDIR}/hf_situation.o

matrix: ${SRCDIR}/matrix.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/matrix.f08 -o ${OUTDIR}/matrix.o

timer: ${SRCDIR}/timer.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/timer.f08 -o ${OUTDIR}/timer.o

mulliken: ${SRCDIR}/mulliken.f08
	${FC} ${FC_FLAGS} -c ${SRCDIR}/mulliken.f08 -o ${OUTDIR}/mulliken.o
