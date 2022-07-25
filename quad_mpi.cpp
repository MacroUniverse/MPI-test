# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <mpi.h>

using namespace std;

void timestamp() {
  const int TIME_SIZE =40;
  static char time_buffer[TIME_SIZE];
  const struct tm *tm_ptr;
  time_t now;
  now = time (NULL);
  tm_ptr = localtime (&now);
  strftime (time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr);
  cout << time_buffer << "\n";
  return;
}

double f(double x) {
  double pi = 3.141592653589793;
  double value = 50.0 / (pi * (2500.0 * x * x + 1.0));
  return value;
}

int main(int argc, char *argv[])
{
  double a = 0, b = 10, error, my_a, my_b, total, wtime, x, my_total;
  int i, n = 10000000, p, q, my_id, my_n, source, tag, target, master = 0;
  MPI_Status status;
  double exact = 0.49936338107645674464;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id); // get process id
  MPI_Comm_size(MPI_COMM_WORLD, &p); // get # of process, only use p-1 to do integral

  cout << "id: " + to_string(my_id) + "  /  " + to_string(p) + "\n";
  if (p == 1) {
	cout << "only 1 process is used! Please run `mpiexec -np 4 ./quad_mpi`. Exiting..." << endl;
	exit(1);
  }

  if (my_id == 0) {
//  If necessary, we adjust N to be divisible by the number of processors.
    my_n = n / (p - 1);
    n = (p - 1) * my_n;

    wtime = MPI_Wtime();

    timestamp();
    cout << "\n";
    cout << "QUAD_MPI\n";
    cout << "  C++/MPI version\n";
    cout << "  Estimate an integral of f(x) from A to B.\n";
    cout << "  f(x) = 50 / (pi * ( 2500 * x * x + 1 ) )\n";
    cout << "\n";
    cout << "  A = " << a << "\n";
    cout << "  B = " << b << "\n";
    cout << "  N = " << n << "\n";
    cout << "  EXACT = " << setw(24) << setprecision(16) << exact << "\n";
    cout << "\n";
    cout << "  Use MPI to divide the computation among\n";
    cout << "  multiple processes.\n";
  }

  source = 0;
  MPI_Bcast(&my_n, 1, MPI_INT, source, MPI_COMM_WORLD);

  // Process 0 assigns each process a subinterval of [A,B].
  if (my_id == 0) {
    for (q = 1; q < p; q++) {
      my_a = ((double)(p - q)*a + (double)(q - 1)*b) / (double)(p - 1);
	  my_b = ((double)(p - q - 1)*a + (double)q*b) / (double)(p - 1);
      target = q;
      tag = 1; int count = 1;
      MPI_Send(&my_a, count, MPI_DOUBLE, target, tag, MPI_COMM_WORLD);
      tag = 2;
      MPI_Send(&my_b, count, MPI_DOUBLE, target, tag, MPI_COMM_WORLD);
    }
    total = 0.0;
    my_total = 0.0;
  }
  else { // Processes receive my_a, my_b, and compute their part of the integral.
    source = 0;
    tag = 1;
    MPI_Recv(&my_a, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
    tag = 2;
    MPI_Recv(&my_b, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
    my_total = 0.0;
    for (i = 1; i <= my_n; i++) {
      x = ( ( double ) ( my_n - i     ) * my_a 
          + ( double ) (        i - 1 ) * my_b )
          / ( double ) ( my_n     - 1 );
      my_total += f(x);
    }
    my_total = (my_b - my_a) / (double)my_n * my_total;
    cout << "  Process " << my_id << " contributed MY_TOTAL = " 
         << my_total << "\n";
  }
//  Each process sends its value to the master process.
  MPI_Reduce(&my_total, &total, 1, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
//  Compute the weighted estimate.
  if (my_id == 0) {
    error = abs(total - exact);
    wtime = MPI_Wtime() - wtime;
    cout << "\n";
    cout << "  Estimate = " << setprecision(16) << total << "\n";
    cout << "  Error = " << error << "\n";
    cout << "  Time = " << wtime << "\n";
  }
//  Terminate MPI.
  MPI_Finalize();
//  Terminate.
  if (my_id == 0) {
    cout << "\n";
    cout << "QUAD_MPI:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp();
  }
  return 0;
}
