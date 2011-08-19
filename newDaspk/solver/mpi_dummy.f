C
C     used only when you cannot access the MPI library or do not want
C     the parallel computation
C
      subroutine mpi_allreduce(
     1     H, HNEW, NUM, MPI_DOUBLE_PRECISION,
     1     MPI_MIN, MPI_COMM_WORLD, IERR)
c=================================================================
c  dummy routine for MPI:
C     collect H from each processor, then choosing the smallest one
C     and broadcast the the new stepsize to each processor
C=================================================================
      DOUBLE PRECISION H, HNEW
      INTEGER MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD
      INTEGER IERR, num
      return
      end
      
