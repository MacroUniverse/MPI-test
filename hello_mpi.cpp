# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>
# include <mpi.h>

using namespace std;

void timestamp()
{
  const int TIME_SIZE = 40;
  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  std::time_t now;
  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );
  std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );
  std::cout << time_buffer << "\n";
  return;
}

// This is a simple MPI test program.
// Each process prints out a "Hello, world!" message.
// The master process also prints out a short message.
// Modified to use the C MPI bindings, 14 June 2016.
int main (int argc, char *argv[])
{
  int id;
  int ierr;
  int p;
  double wtime;
  ierr = MPI_Init(&argc, &argv);

  if (ierr != 0) {
    cout << "\n";
    cout << "HELLO_MPI - Fatal error!\n";
    cout << "  MPI_Init returned nonzero ierr.\n";
    exit(1);
  }
  // Get the number of processes.
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &p);
  if (p == 1) {
    cout << "only 1 process is used! Please run `mpiexec -np 4 ./hello_mpi`. Exiting..." << endl;
    exit(1);
  }
  // Get the individual process ID.
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &id);
  // Process 0 prints an introductory message.
  if (id == 0) {
    timestamp ( );
    cout << "\n";
    cout << "P" << id << ":  HELLO_MPI - Master process:\n";
    cout << "P" << id << ":    C++/MPI version\n";
    cout << "P" << id << ":    An MPI example program.\n";
    cout << "\n";
    cout << "P" << id << ":    The number of processes is " << p << "\n";
    cout << "\n";
  }
//  Every process prints a hello.
  if (id == 0)
    wtime = MPI_Wtime();
  cout << "P" << id << ":    'Hello, world!'\n";
//  Process 0 says goodbye.
  if (id == 0) {
    wtime = MPI_Wtime() - wtime;
    cout << "P" << id << ":    Elapsed wall clock time = " << wtime << " seconds.\n";
  }
  MPI_Finalize();
//  Terminate.
  if (id == 0) {
    cout << "\n";
    cout << "P" << id << ":  HELLO_MPI:\n";
    cout << "P" << id << ":  Normal end of execution.\n";
    cout << "\n";
    timestamp();
  }
  return 0;
}
