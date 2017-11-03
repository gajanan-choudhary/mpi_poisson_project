#ifdef _MESSG
#ifdef _MPI
  void message_error(int);
  void message_barrier(MPI_Comm);
#endif
#endif
