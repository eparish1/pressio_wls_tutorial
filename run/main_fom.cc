
#include "../swe_app/swe_app.hpp"
void solve_FOM(double g, double mu, int sol_indx){
  int nx = 128;   
  int ny = 128;  
  int N_cell = nx*ny; 
  double L = 5.;
  double dx = L/(nx );
  double dy = L/(ny );
  double dt = 0.005;
  double et = 5.;
  double t = 0;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
    using fom_t = ShallowWaterApp<double>;
    using scalar_t        = typename fom_t::scalar_type;
    using native_state_t  = typename fom_t::state_type;
    using native_dmat_t      = typename fom_t::dense_matrix_type;
    using fom_state_t     = pressio::containers::Vector<native_state_t>;
    using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;

    using execution_space = Kokkos::DefaultExecutionSpace;
    using kll		= Kokkos::LayoutLeft;
    using k1dLl_d		= Kokkos::View<scalar_t*, kll, execution_space>;
    using k2dLl_d		= Kokkos::View<scalar_t**, kll, execution_space>;
    using rom_state_t	= pressio::containers::Vector<k1dLl_d>;
    using hessian_t	= pressio::containers::Matrix<k2dLl_d>;
    using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;
    using ode_tag   = ::pressio::ode::implicitmethods::Euler;

    fom_t appObj(nx,ny,dx,dy,romSize,g,Comm);
    fom_t::state_t U = appObj.getGaussianIC(mu);
    fom_t::state_t V = appObj.createVelocity();
    fom_t::state_t U0 = appObj.createVelocity();

    double rk4const[ 4 ];
    rk4const[0] = 1./4.;
    rk4const[1] = 2./4.;
    rk4const[2] = 3./4.;
    rk4const[3] = 1.;
    int save_freq = 10;
    int counter = 0;
    std::string filename = "solution";
    filename += std::to_string(sol_indx);
    filename += ".bin";
    std::ofstream myfile (filename,  std::ios::out | std::ios::binary); 

    while (t <= et - dt/2.){
      fom_t::state_t U0(U,Teuchos::DataAccess::Copy);
      if (counter%save_freq == 0){
        for (int i = 0;i < N_cell; i++){
          double *uLocal;
          U.getLocalRowView(i,uLocal);
          auto v = uLocal[0];
          myfile.write(reinterpret_cast<const char*>(&uLocal[0]),sizeof(uLocal[0]));
          myfile.write(reinterpret_cast<const char*>(&uLocal[1]),sizeof(uLocal[1]));
          myfile.write(reinterpret_cast<const char*>(&uLocal[2]),sizeof(uLocal[2]));

        }
      }
   
      for (int k = 0; k < 4; k++){
        appObj.velocity(U,0.,V);
        for (int i=0; i < N_cell; i++){
          double *Up;
          double *U0p;
          double *Vp;
          U.getLocalRowView(i,Up);
          U0.getLocalRowView(i,U0p);
          V.getLocalRowView(i,Vp);
          for (int j=0;j<3;j++){
            Up[j] = U0p[j] + dt*rk4const[k]*Vp[j];
          } 
        }
      }
      t += dt;
      counter += 1;
      std::cout << t << std::endl;
    }
    myfile.close(); 
}



int main( int argc, char* argv[] )
{
  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  std::string checkStr {"PASSED"};
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {

    std::array<double,3> ga = {3.,6.,9.};
    std::array<double,4> mua = {0.05,0.1,0.2};
    int counter = 0;
    
    for (int i=0; i < 3; i++){
      for (int j=0; j < 3; j++){
        auto startTime = std::chrono::high_resolution_clock::now();
        solve_FOM(ga[i],mua[j],counter);  
        counter += 1; 
        const auto finishTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> elapsed = finishTime - startTime;
        std::cout << "Walltime = " << elapsed.count() << '\n';
      }
    } 
    
    double g_test = 7.5;
    double mu_test = 0.125;
    solve_FOM(g_test,mu_test,100);
  }
  return 0;
}
//}}}

