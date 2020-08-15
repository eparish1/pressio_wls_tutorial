#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <iostream>
#include <fstream>

template<typename flux_t,typename state_t, typename normal_t>
void roeflux_kernel( flux_t  & F, const state_t & UL, const state_t & UR, const normal_t n, double g){
// PURPOSE: This function calculates the flux for the Euler equations
// using the Roe flux function
// INPUTS:
//    UL: conservative state vector in left cell
//   UR: conservative state vector in right cell
//    n: normal pointing from the left cell to the right cell

// OUTPUTS:
//  F   : the flux out of the left cell (into the right cell)
  double es = 1.e-30;
  double hL = UL[0];
  double uL = UL[1]/(hL + es);
  double vL = UL[2]/(hL + es);
  double unL = uL*n[0] + vL*n[1];

  double pL = 0.5*g*pow(hL,2.);
  double FL[ 3 ];
  FL[0] = hL*unL;
  FL[1] = UL[1]*unL + pL*n[0];
  FL[2] = UL[2]*unL + pL*n[1];

  double hR = UR[0];
  double uR = UR[1]/(hR + es);
  double vR = UR[2]/(hR + es);
  double unR = uR*n[0] + vR*n[1];
  double pR = 0.5*g*pow(hR,2.);
  double FR[3];
  FR[0] = hR*unR;
  FR[1] = UR[1]*unR + pR*n[0];
  FR[2] = UR[2]*unR + pR*n[1];

  // rho average
  double hm = 0.5*(hL + hR);
  double um = (unL*pow(hL,0.5) + unR*pow(hR,0.5) )/( pow(hL,0.5) + pow(hR,0.5)  + es);
  // eigenvalues
  double smax = abs(um) + abs(pow(g*hm, 0.5));
  // flux assembly
  F[0]    = 0.5*(FL[0]+FR[0])-0.5*smax*(UR[0] - UL[0]);
  F[1]    = 0.5*(FL[1]+FR[1])-0.5*smax*(UR[1] - UL[1]);
  F[2]    = 0.5*(FL[2]+FR[2])-0.5*smax*(UR[2] - UL[2]);

}



template<typename scalar_t>
class ShallowWaterAppHyper
{
  protected:
    using map_t		= Tpetra::Map<>;
    using nativeVec	= Tpetra::BlockVector<>;
 
    using go_t		= typename map_t::global_ordinal_type;
    using lo_t		= typename map_t::local_ordinal_type;
  
    using tcomm_t		= Teuchos::MpiComm<int>;
    using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
    using rcpmap_t	= Teuchos::RCP<const map_t>;

    template<typename T> using stdrcp = std::shared_ptr<T>;
    using crs_graph_type = Tpetra::CrsGraph<>;

  private:
    scalar_t dx_;
    scalar_t dy_;
    int nx_;
    int ny_;
    scalar_t g_;


  protected:
    Teuchos::RCP<Teuchos::FancyOStream> wrappedCout_;
    rcpcomm_t comm_{};
    rcpmap_t hyperMap_{};
    rcpmap_t hyperMapWithStencil_{};

    int myRank_{};
    int totRanks_{};

    mutable stdrcp<nativeVec> U_{}; 
    stdrcp<Tpetra::Vector<>> xGrid_{}; 
    stdrcp<Tpetra::Vector<>> yGrid_{};

    std::vector<go_t> sm_gids_; 
    std::vector<go_t> smps_gids_; 

    int sampleMeshSize_{};
    int sampleMeshPlusStencilSize_{};

  public:
    using scalar_type	= scalar_t;
    using state_type	= nativeVec;
    using state_t = state_type;
    using velocity_type	= state_type;
    using jacobian_type	= Tpetra::CrsMatrix<>;
    using dense_matrix_type = Tpetra::BlockMultiVector<>;

    
  public:
    ShallowWaterAppHyper(int nx, int ny, scalar_t dx, scalar_t dy, int sampleMeshSize, int sampleMeshPlusStencilSize, double g, rcpcomm_t comm) : 
      nx_(nx),dx_(dx),ny_(ny),dy_(dy), sampleMeshSize_(sampleMeshSize), sm_gids_(sampleMeshSize), sampleMeshPlusStencilSize_(sampleMeshPlusStencilSize),
      smps_gids_(sampleMeshPlusStencilSize), g_(g), comm_(comm){this->setup();}

    std::array<int,2> get_ij_from_gid(int gid) const {
      int j = gid/nx_;
      int i = gid%nx_;
      std::array<int,2> ij = {i,j};
      return ij;
    }

    int modulus(int n,int M) const {
      return ((n % M) + M) % M;
    }

    int get_gid_from_ij(int i,int j) const {
      int global_indx = (modulus(j,ny_))*nx_ + modulus(i,nx_);
      return global_indx;
    } 

    //========== 
    void velocity(const state_t & U, const scalar_t & t, velocity_type & V) const {

      double *UL;
      double *UR;
      double *UU;
      double *UD;
      double *V_view;

      double FL[3];
      double FR[3];
      double FU[3];
      double FD[3];

      std::array<double, 2> nx = { 1, 0 };
      std::array<double, 2> ny = { 0, 1 };

      const auto gids = hyperMap_->getMyGlobalIndices();
      //================= 
      for (int sid=0; sid < sampleMeshSize_  ; sid++){
        auto gid = gids[sid];
        auto ij = get_ij_from_gid(gid);
        U.getGlobalRowView(get_gid_from_ij(ij[0] - 1,ij[1]),UL);
        U.getGlobalRowView(get_gid_from_ij(ij[0]    ,ij[1]),UR);
        roeflux_kernel(FL,UL,UR,nx,g_);

        U.getGlobalRowView(get_gid_from_ij(ij[0]    , ij[1]),UL);
        U.getGlobalRowView(get_gid_from_ij(ij[0] + 1, ij[1]),UR);
        roeflux_kernel(FR,UL,UR,nx,g_);

        U.getGlobalRowView(get_gid_from_ij(ij[0]  ,ij[1] - 1),UD);
        U.getGlobalRowView(get_gid_from_ij(ij[0]  ,ij[1]    ),UU);
        roeflux_kernel(FD,UD,UU,ny,g_);

        U.getGlobalRowView(get_gid_from_ij(ij[0]  ,ij[1]    ),UD);
        U.getGlobalRowView(get_gid_from_ij(ij[0]  ,ij[1] + 1),UU);
        roeflux_kernel(FU,UD,UU,ny,g_);


        V.getGlobalRowView(gid,V_view); 
        for (int j=0;j<3;j++){
          V_view[j] = -1./dx_*(FR[j] - FL[j]) - 1./dy_*(FU[j] - FD[j]);
        }
      }
    }

    
    void applyJacobian(const state_t &U, const dense_matrix_type & A,scalar_t t, dense_matrix_type &JA)const{
      double eps = 1.e-5;
      state_t Up(U,Teuchos::DataAccess::Copy);
      velocity_type V0(*hyperMap_,3);
      velocity_type V_perturb(*hyperMap_,3);
      V_perturb.putScalar(0.);

      velocity(U,t,V0);
      const auto gids = hyperMap_->getMyGlobalIndices();

      for (int k=0; k < A.getNumVectors(); k++){

        for (int i=0; i < sampleMeshSize_  ; i++){

          double *Uview = nullptr;
          double *Upview = nullptr;
          double *Aview = nullptr;
          auto gid = gids[i];
          U.getGlobalRowView(gid,Uview);
          Up.getGlobalRowView(gid,Upview);
          A.getGlobalRowView(gid,k,Aview);
          for (int j=0;j<3;j++){
            Upview[j] = Uview[j] + eps*Aview[j];
          }
        }
        velocity(Up,t,V_perturb);

        for (int i=0; i < sampleMeshSize_  ; i++){
          double *V0_view = nullptr;
          double *V_perturb_view = nullptr;
          double *JA_view = nullptr;
          auto gid = gids[i];
          V0.getGlobalRowView(gid,V0_view);
          V_perturb.getGlobalRowView(gid,V_perturb_view);
          JA.getGlobalRowView(gid,k,JA_view);
          for (int j=0; j < 3; j++){
            JA_view[j] = 1./eps*(V_perturb_view[j] - V0_view[j] ) ;
          }
        }

      }
    }


    //========
    velocity_type createVelocity() const {
      velocity_type V(*hyperMap_,3);
      return V;
    }

    dense_matrix_type createApplyJacobianResult(const dense_matrix_type & A) const
    {
      dense_matrix_type JA(*hyperMap_, 3, A.getNumVectors() );
      return JA;
    }


    rcpmap_t getHyperMapWithStencil(){
      return hyperMapWithStencil_;
    };

protected:
  void setup(){
    using Teuchos::rcpFromRef;
    using Teuchos::FancyOStream;
    wrappedCout_ = getFancyOStream(rcpFromRef(std::cout));

    myRank_ =  comm_->getRank();
    totRanks_ =  comm_->getSize();

    std::ifstream sample_mesh_gids_file("sample_mesh_gids.txt");
 
    for (int i =0; i < sampleMeshSize_ ; i++){
      sample_mesh_gids_file >> sm_gids_[i];
    }
    sample_mesh_gids_file.close();

    std::ifstream sample_mesh_gids_plus_stencil_file("sample_mesh_plus_stencil_gids.txt");
    for (int i =0; i < sampleMeshPlusStencilSize_; i++){
      sample_mesh_gids_plus_stencil_file >>smps_gids_[i];
    }

    sample_mesh_gids_plus_stencil_file.close(); 
  
    hyperMap_ = Tpetra::createNonContigMap<lo_t,go_t>(sm_gids_,comm_);
    hyperMapWithStencil_ = Tpetra::createNonContigMap<lo_t,go_t>(smps_gids_,comm_);

    U_ = std::make_shared<nativeVec>(*hyperMapWithStencil_,3);

    for (int i=0; i < sampleMeshSize_; i++){
      auto id = hyperMap_->getGlobalElement(i);   
    }


    xGrid_ = std::make_shared<Tpetra::Vector<>>(hyperMapWithStencil_);
    yGrid_ = std::make_shared<Tpetra::Vector<>>(hyperMapWithStencil_);
    for (int sid = 0; sid < sampleMeshPlusStencilSize_; sid++){
        auto gid = hyperMapWithStencil_->getGlobalElement(sid); 
        auto ij = get_ij_from_gid(gid);
        xGrid_->replaceGlobalValue(gid, dx_*ij[0] + dx_*0.5);
        yGrid_->replaceGlobalValue(gid, dy_*ij[1] + dy_*0.5);
    }

    };

public:
 
  nativeVec getGaussianIC(double mu){
    auto xGridv = xGrid_->getData();
    auto yGridv = yGrid_->getData();
    int i = 0;
    for (int i=0; i < sampleMeshPlusStencilSize_; i++){
        double * uVal;
        U_->getLocalRowView(i,uVal);
        uVal[0] = 1. + mu*exp( - ( pow( xGridv[i] - 2.5,2)  + pow(yGridv[i] - 2.5,2) ));
        uVal[1] = 0.;
        uVal[2] = 0.; 
    }
    return *U_;
  }

};





template<typename scalar_t>
class ShallowWaterApp
{
  protected:
    using map_t		= Tpetra::Map<>;
    using nativeVec	= Tpetra::BlockVector<>;
 
    using go_t		= typename map_t::global_ordinal_type;
    using lo_t		= typename map_t::local_ordinal_type;
  
    using tcomm_t		= Teuchos::MpiComm<int>;
    using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
    using rcpmap_t	= Teuchos::RCP<const map_t>;

    template<typename T> using stdrcp = std::shared_ptr<T>;
    using crs_graph_type = Tpetra::CrsGraph<>;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> wrappedCout_;
    rcpcomm_t comm_{};
    rcpmap_t gridMap_{};
 
    int myRank_{};
    int totRanks_{};
    lo_t NumMyElem_{};
    std::vector<go_t> myGel_{};
  
    mutable stdrcp<nativeVec> U_{}; 
    mutable stdrcp<nativeVec> U_tmp_{}; 
    mutable stdrcp<nativeVec> V_tmp_{}; 
    mutable stdrcp<nativeVec> V_tmp2_{}; 
    stdrcp<Tpetra::Vector<>> xGrid_{}; // mesh points coordinates
    stdrcp<Tpetra::Vector<>> yGrid_{}; // mesh points coordinates


  private:
    scalar_t dx_;
    scalar_t dy_;

    int nx_;
    int ny_;
    double g_;
  public:
    using scalar_type	= scalar_t;
    using state_type	= nativeVec;
    using state_t = state_type;
    using velocity_type	= state_type;
    using jacobian_type	= Tpetra::CrsMatrix<>;
    using dense_matrix_type = Tpetra::BlockMultiVector<>;

    
  public:
    ShallowWaterApp(int nx, int ny, scalar_t dx, scalar_t dy, double g, rcpcomm_t comm) : 
      nx_(nx),dx_(dx),ny_(ny),dy_(dy),g_(g), comm_(comm){this->setup();}


    int modulus(int n,int M) const {
      return ((n % M) + M) % M;
    }

    int index_mapper(int i,int j) const {
 
      int global_indx = (modulus(j,ny_))*nx_ + modulus(i,nx_);
      return global_indx;
    } 
 
    //========== 
    void velocity(const state_t & U, const scalar_t & t, velocity_type & V) const {
      double *UL;
      double *UR;
      double *UU;
      double *UD;
      double *V_view;

      double FL[3];
      double FR[3];
      double FU[3];
      double FD[3];

      std::array<double, 2> nx = { 1, 0 };
      std::array<double, 2> ny = { 0, 1 };

      
      for (int i=0; i < nx_  ; i++){
        for (int j=0; j < ny_ ; j++){

          U.getLocalRowView(index_mapper(i-1,j),UL);
          U.getLocalRowView(index_mapper(i,j),UR);
          auto indx= index_mapper(i-1,j);
          roeflux_kernel(FL,UL,UR,nx,g_);
          U.getLocalRowView(index_mapper(i,j),UL);
          U.getLocalRowView(index_mapper(i+1,j),UR);
          roeflux_kernel(FR,UL,UR,nx,g_);

          U.getLocalRowView(index_mapper(i,j-1),UD);
          U.getLocalRowView(index_mapper(i,j),UU);
          roeflux_kernel(FD,UD,UU,ny,g_);
          U.getLocalRowView(index_mapper(i,j),UD);
          U.getLocalRowView(index_mapper(i,j+1),UU);
          roeflux_kernel(FU,UD,UU,ny,g_);

          V.getLocalRowView(index_mapper(i,j),V_view); 
          for (int k=0;k<3;k++){
            V_view[k] = -1./dx_*(FR[k] - FL[k]) - 1./dy_*(FU[k] - FD[k]);
          }
        }
      }
    }
    /*

    */
    
    void applyJacobian(const state_t &U, const dense_matrix_type & A,scalar_t t, dense_matrix_type &JA)const{
      double eps = 1.e-5;
      state_t Up(*gridMap_,3);
      velocity_type V0(*gridMap_,3);
      velocity_type V_perturb(*gridMap_,3);

      velocity(U,t,V0);
      
      for (int k=0; k < A.getNumVectors(); k++){
        int i = 0;
        for (auto const & it : myGel_){
          double *Uview;
          double *Upview;
          double *Aview;
          U.getLocalRowView(i,Uview);
          Up.getLocalRowView(i,Upview);
          A.getLocalRowView(i,k,Aview);
          for (int j=0;j<3;j++){
            Upview[j] = Uview[j] + eps*Aview[j];
          }
          i++;
        }
        velocity(Up,t,V_perturb);
        i = 0;
        for (auto const & it : myGel_){
            double *V0_view;
            double *V_perturb_view;
            double *JA_view;
            V0.getLocalRowView(i,V0_view);
            V_perturb.getLocalRowView(i,V_perturb_view);
            JA.getLocalRowView(i,k,JA_view);
            for (int j=0; j < 3; j++){
              JA_view[j] = 1./eps*(V_perturb_view[j] - V0_view[j] ) ;
            }
          i++;
          }
        }
      }

    //========
    velocity_type createVelocity() const {
      velocity_type V(*gridMap_,3);
      return V;
    }

    dense_matrix_type createApplyJacobianResult(const dense_matrix_type & A) const
    {
      dense_matrix_type JA(*gridMap_, 3, A.getNumVectors() );
      return JA;
    }


    rcpmap_t getGridMap(){
      return gridMap_;
    };


   int getNumMyElem() const{
     return NumMyElem_;
   }
protected:
  void setup(){
    using Teuchos::rcpFromRef;
    using Teuchos::FancyOStream;
    wrappedCout_ = getFancyOStream(rcpFromRef(std::cout));

    myRank_ =  comm_->getRank();
    totRanks_ =  comm_->getSize();

    // distribute cells
    int N_cell = nx_*ny_;
    gridMap_ = Teuchos::rcp(new map_t(N_cell, 0, comm_));

    xGrid_ = std::make_shared<Tpetra::Vector<>>(gridMap_);
    yGrid_ = std::make_shared<Tpetra::Vector<>>(gridMap_);
    auto xGridv = xGrid_->getDataNonConst();
    auto yGridv = yGrid_->getDataNonConst();

    for (int i=0; i < nx_; i ++){
      for (int j=0; j < ny_; j++){
        xGridv[index_mapper(i,j)] = dx_*i + dx_*0.5;
        yGridv[index_mapper(i,j)] = dy_*j + dy_*0.5;
      }
    }
    U_ = std::make_shared<nativeVec>(*gridMap_,3);
    U_->putScalar(0.0);

    U_tmp_ = std::make_shared<nativeVec>(*gridMap_,3);
    U_tmp_->putScalar(0.0);

    V_tmp_ = std::make_shared<nativeVec>(*gridMap_,3);
    V_tmp_->putScalar(0.0);

    V_tmp2_= std::make_shared<nativeVec>(*gridMap_,3);
    V_tmp2_->putScalar(0.0);

  };

public:
  nativeVec getGaussianIC(double mu){
    auto xGridv = xGrid_->getData();
    auto yGridv = yGrid_->getData();
    int i = 0;
    for (int i=0; i < nx_; i++){
      for (int j=0; j < ny_; j++){
        double * uVal;
        U_->getLocalRowView(index_mapper(i,j),uVal);
        uVal[0] = 1. + mu*exp( - ( pow( xGridv[index_mapper(i,j)] - 2.5,2)  + pow(yGridv[index_mapper(i,j)] - 2.5,2) ));
        uVal[1] = 0.;
        uVal[2] = 0.; 
      }
    }
    return *U_;
  }

};


