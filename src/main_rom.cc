
#include "pressio_rom.hpp"
#include "swe_hyper.hpp"

void solve_ROM(double g, double mu, int sol_indx)
{
  int nx = 128;
  int ny = 128;
  int N_cell = nx*ny;
  double L = 5.;
  double dx = L/(nx );
  double dy = L/(ny );
  double dt = 0.05;
  double et = 5.;
  double t = 0;

  using tcomm_t		= Teuchos::SerialComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;

  // int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
  rcpcomm_t Comm = Teuchos::rcp (new tcomm_t());
  using fom_t = ShallowWaterAppHyper<double>;
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
  using ode_tag   = ::pressio::ode::implicitmethods::BDF2;

  int romSize;
  int sampleMeshSize;
  int sampleMeshPlusStencilSize;
  std::ifstream info_file("info_file.txt");
  info_file >> romSize;
  std::cout << romSize << std::endl;;
  info_file >> sampleMeshSize;
  info_file >> sampleMeshPlusStencilSize;
  info_file.close();
  std::cout << "| romSize = " << romSize
	    << "| |  sampleMeshSize = " << sampleMeshSize
	    << "| |  sampleMeshPlusStencilSize " << sampleMeshPlusStencilSize
	    << std::endl;
  
  fom_t appObj(nx,ny,dx,dy,sampleMeshSize,sampleMeshPlusStencilSize,g,Comm);
  fom_t::state_t U = appObj.getGaussianIC(mu);
  fom_t::state_t V = appObj.createVelocity();

  std::ifstream file("PhiSamplePlusStencil.txt");
  auto hyperMapWithStencil = appObj.getHyperMapWithStencil();
  fom_t::dense_matrix_type Phi(*hyperMapWithStencil,3,romSize);
  for (int i =0; i < sampleMeshPlusStencilSize ; i++){
    for (int j=0 ; j < 3 ; j++){
      for (int k=0; k < romSize ; k++){
	double * PhiView;
	auto gid = hyperMapWithStencil->getGlobalElement(i);
	Phi.getGlobalRowView(gid,k,PhiView);
	file >> PhiView[j];
      }
    }
  }

  file.close();

  const fom_state_t fomStateInitCond(U);
  fom_state_t fomStateReference(U);
  pressio::ops::set_zero(fomStateReference);
  pressio::rom::wls::window_size_t numStepsInWindow = 2;
  pressio::rom::wls::rom_size_t wlsSize = romSize*numStepsInWindow;
  scalar_t finalTime = 5;
  pressio::rom::wls::window_size_t numWindows = (finalTime/dt)/numStepsInWindow;
  decoder_t decoderObj(Phi);
  
  //  lin solver
  using lin_solver_tag  = pressio::solvers::linear::direct::potrsL;
  using linear_solver_t = pressio::solvers::linear::Solver<lin_solver_tag, hessian_t>;
  linear_solver_t linear_solver;

  //WLS problem ***
  using precon_type = ::pressio::rom::wls::preconditioners::NoPreconditioner;
  using jacobians_update_tag = ::pressio::rom::wls::FrozenJacobian;
  using policy_t     = pressio::rom::wls::HessianGradientSequentialPolicy<
    fom_t, decoder_t,ode_tag, pressio::matrixLowerTriangular,
    precon_type, jacobians_update_tag>;

  using wls_system_t = pressio::rom::wls::SystemHessianAndGradientApi<
    rom_state_t, decoder_t, ode_tag, hessian_t, policy_t>;

  // create policy and wls system
  int jacobianUpdateFrequency = 10;
  policy_t hgPolicy(romSize, numStepsInWindow, decoderObj,
		    appObj, fomStateReference,
		    wls_system_t::timeStencilSize_,jacobianUpdateFrequency);
  wls_system_t wlsSystem(romSize, numStepsInWindow,
			 decoderObj, hgPolicy, fomStateInitCond,
			 fomStateReference, linear_solver);
  // create the wls state
  rom_state_t  wlsState(wlsSize);

  using gn_t = pressio::solvers::nonlinear::composeGaussNewton_t<
    wls_system_t,
    pressio::solvers::nonlinear::DefaultUpdate,
    linear_solver_t>;

  gn_t GNSolver(wlsSystem, wlsState, linear_solver);
  GNSolver.setTolerance(1e-3);
  GNSolver.setMaxIterations(100);

  //  solve wls problem
  auto startTime = std::chrono::high_resolution_clock::now();
  fom_state_t fomStatePP(U);
  std::string filename = "wls_rom_solution";
  filename += std::to_string(sol_indx);
  filename += ".txt";

  std::ofstream myfile (filename);
  int save_freq = 1;
  int counter = 0;
  for (auto iWind = 0; iWind < numWindows; iWind++){
    wlsSystem.advanceOneWindow(wlsState, GNSolver, iWind, dt);
    if (counter%save_freq == 0){
      auto wlsData = *wlsState.data();
      for (int i = 0; i < romSize*numStepsInWindow; i++){
	auto indx = i;
	myfile << wlsData(indx) << std::endl;
      }
    }
    counter += 1*numStepsInWindow;
  }
  myfile.close();
  
}


int main( int argc, char* argv[] )
{
  using tcomm_t		= Teuchos::SerialComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  std::string checkStr {"PASSED"};
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    auto startTime = std::chrono::high_resolution_clock::now();
    double g_test = 7.5;
    double mu_test = 0.125;
    solve_ROM(g_test,mu_test,100);
    const auto finishTime = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = finishTime - startTime;
    std::cout << "Walltime = " << elapsed.count() << '\n';

  }
  return 0;
}
