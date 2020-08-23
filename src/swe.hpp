
#ifndef SWE_HPP_
#define SWE_HPP_

#include "roeflux.hpp"
#include <Tpetra_BlockVector_decl.hpp>
#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_Map_decl.hpp>

template<typename scalar_t>
class ShallowWaterApp
{

protected:
  using map_t		= Tpetra::Map<>;
  using nativeVec	= Tpetra::BlockVector<>;

  using go_t		= typename map_t::global_ordinal_type;
  using lo_t		= typename map_t::local_ordinal_type;

  using tcomm_t		= Teuchos::SerialComm<int>;
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
  ShallowWaterApp(int nx, int ny,
		  scalar_t dx, scalar_t dy, double g,
		  rcpcomm_t comm) :
    nx_(nx),dx_(dx),ny_(ny),dy_(dy),g_(g), comm_(comm){this->setup();}


  int modulus(int n,int M) const {
    return ((n % M) + M) % M;
  }

  int index_mapper(int i,int j) const {

    int global_indx = (modulus(j,ny_))*nx_ + modulus(i,nx_);
    return global_indx;
  }

  //==========
  void velocity(const state_t & U, const scalar_t & t, velocity_type & V) const
  {
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

  void applyJacobian(const state_t &U,
		     const dense_matrix_type & A,
		     scalar_t t,
		     dense_matrix_type &JA)const
  {
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
        uVal[0] = 1. + mu*exp( - ( pow( xGridv[index_mapper(i,j)] - 2.5,2)
				   + pow(yGridv[index_mapper(i,j)] - 2.5,2) ));
        uVal[1] = 0.;
        uVal[2] = 0.;
      }
    }
    return *U_;
  }

};

#endif
