
#ifndef SWE_HYPER_HPP_
#define SWE_HYPER_HPP_

#include "roeflux.hpp"
#include <Tpetra_BlockMultiVector_decl.hpp>
#include <Tpetra_BlockMultiVector_def.hpp>
#include <Tpetra_BlockMultiVector_decl.hpp>
#include <Tpetra_BlockMultiVector_fwd.hpp>
#include "/Users/ejparis/tutorial_build/tpls/trilinos/Trilinos/packages/tpetra/core/src/Tpetra_BlockMultiVector_decl.hpp"
//#include 
#include <Tpetra_BlockVector_decl.hpp>
#include <Tpetra_Experimental_BlockMultiVector_decl.hpp>
#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>

template<typename scalar_t>
class ShallowWaterAppHyper
{
  protected:
    using map_t		= Tpetra::Map<>;
    using nativeVec	= Tpetra::BlockVector<>;

    using go_t		= typename map_t::global_ordinal_type;
    using lo_t		= typename map_t::local_ordinal_type;

    using tcomm_t	= Teuchos::SerialComm<int>;
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

  public:
    using scalar_type	= scalar_t;
    using state_type	= nativeVec;
    using state_t = state_type;
    using velocity_type	= state_type;
    using jacobian_type = Tpetra::BlockCrsMatrix<>;
    using dense_matrix_type = Tpetra::BlockMultiVector<>;



  protected:
    Teuchos::RCP<Teuchos::FancyOStream> wrappedCout_;
    rcpcomm_t comm_{};
    rcpmap_t hyperMap_{};
    rcpmap_t hyperMapWithStencil_{};
    rcpmap_t hyperMapWithStencil2_{};
    using importer_t = Tpetra::Import<lo_t,go_t>;
    stdrcp<importer_t>	importer_{};
    int myRank_{};
    int totRanks_{};

    mutable stdrcp<nativeVec> U_{};

    stdrcp<jacobian_type> Jac_{};
    stdrcp<Tpetra::Vector<>> xGrid_{};
    stdrcp<Tpetra::Vector<>> yGrid_{};

    std::vector<go_t> sm_gids_;
    std::vector<go_t> smps_gids_;

    int sampleMeshSize_{};
    int sampleMeshPlusStencilSize_{};



  public:
    ShallowWaterAppHyper(int nx, int ny,
			 scalar_t dx, scalar_t dy,
			 int sampleMeshSize, int sampleMeshPlusStencilSize,
			 double g, rcpcomm_t comm)
      : nx_(nx),dx_(dx),ny_(ny),dy_(dy),
	sampleMeshSize_(sampleMeshSize),
	sm_gids_(sampleMeshSize),
	sampleMeshPlusStencilSize_(sampleMeshPlusStencilSize),
	smps_gids_(sampleMeshPlusStencilSize),
	g_(g),
	comm_(comm){
      this->setup();
    }

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


    
    //void updateJacobian(const state_t & U, const scalar_t & t, velocity_type & V) const {
    void updateJacobian(const state_type & U,
  			     jacobian_type & jac,
			     const scalar_type ) const{

      // J is of size sampleMeshSize_ x sampleMeshSize_
      double *UL;
      double *UR;
      double *UU;
      double *UD;
      double *V_view;

      double JL_L[3][3];
      double JR_L[3][3];
      double JL_R[3][3];
      double JR_R[3][3];
      double JD_D[3][3];
      double JU_D[3][3];
      double JD_U[3][3];
      double JU_U[3][3];

      std::array<double, 2> nx = { 1, 0 };
      std::array<double, 2> ny = { 0, 1 };


      const auto gids = hyperMap_->getMyGlobalIndices();
      //=================
      for (int sid=0; sid < sampleMeshSize_  ; sid++){
        auto gid = gids[sid];
        auto ij = get_ij_from_gid(gid);
        auto gid_im1 = get_gid_from_ij(ij[0] - 1,ij[1]);
        auto gid_ip1 = get_gid_from_ij(ij[0] + 1,ij[1]);
        auto gid_jm1 = get_gid_from_ij(ij[0],ij[1] - 1);
        auto gid_jp1 = get_gid_from_ij(ij[0],ij[1] + 1);

        U.getGlobalRowView(gid_im1,UL);
        U.getGlobalRowView(gid,UR);
        roeflux_jacobian(JL_L,JR_L,UL,UR,nx,g_);

        U.getGlobalRowView(gid,UL);
        U.getGlobalRowView(gid_ip1,UR);
        roeflux_jacobian(JL_R,JR_R,UL,UR,nx,g_);

        U.getGlobalRowView(gid_jm1,UD);
        U.getGlobalRowView(gid,UU);
        roeflux_jacobian(JD_D,JU_D,UD,UU,ny,g_);

        U.getGlobalRowView(gid,UD);
        U.getGlobalRowView(gid_jp1,UU);
        roeflux_jacobian(JD_U,JU_U,UD,UU,ny,g_);

        // column index is gid
        // row index
        lo_t blockSize = 3;
        const lo_t * colInds;
        scalar_t* vals;

    	  int err = 0;
        int nvals = jac.getNumEntriesInLocalRow(sid);
        std::vector<go_t> col_gids(nvals);
        err = jac.getLocalRowView(sid,colInds, vals,nvals);
        if (err != 0) {
         std::cout << "error in getting row view" << std::endl;
         break;
        }
        auto colMap = jac.getColMap();
        for (lo_t col_no = 0; col_no < nvals; col_no ++){
          col_gids[col_no] = colMap->getGlobalElement(colInds[col_no]); 
          //std::cout << sid << " " << col_no << " " << col_gids[col_no] << std::endl;
        }
       /// 
       for (int k = 0; k < nvals; k++){
          auto col_gid = col_gids[k];

      if (col_gid == gid){
        for (lo_t j = 0; j < blockSize; ++j) {
          for (lo_t i = 0; i < blockSize; ++i) {
            vals[blockSize * blockSize * k + i + j * blockSize] = -1./dx_*(JL_R[i][j] - JR_L[i][j]) - 1./dy_*(JD_U[i][j] - JU_D[i][j]) ;
          }
        }
      }
      else if (col_gid == gid_im1){
        // Blocks are stored in row-major format.
        for (lo_t j = 0; j < blockSize; ++j) {
          for (lo_t i = 0; i < blockSize; ++i) {
            vals[blockSize * blockSize * k + i + j * blockSize] = 1./dx_*JL_L[i][j] ;
          }
        }
      }
      else if (col_gid == gid_ip1){
        for (int j = 0; j < blockSize; ++j) {
          for (int i = 0; i < blockSize; ++i) {
            vals[blockSize * blockSize * k + i + j * blockSize] = -1./dx_*JR_R[i][j] ;
          }
        }
      }

      else if (col_gid == gid_jm1){
        for (int j = 0; j < blockSize; ++j) {
          for (int i = 0; i < blockSize; ++i) {
            vals[blockSize * blockSize * k + i + j * blockSize] = 1./dy_*JD_D[i][j] ;
          }
        }
      }

      else if (col_gid == gid_jp1){
        for (int j = 0; j < blockSize; ++j) {
          for (int i = 0; i < blockSize; ++i) {
            vals[blockSize * blockSize * k + i + j * blockSize] = -1./dy_*JU_U[i][j] ;
          }
        }
      }
        }

      }

    } 


    // computes: A = Jac B where B is dense
    void applyJacobian(const state_type & y,
           const dense_matrix_type & B, //of size sample mesh plus stencil
           scalar_type t,
           dense_matrix_type & A) const //of size sample mesh
    {
      dense_matrix_type B_sm(*hyperMap_, 3, B.getNumVectors() );
      auto  B_sm_vv = B_sm.getMultiVectorView();
      auto  B_vv = B.getMultiVectorView();
      B_sm_vv.doImport(B_vv,*importer_,Tpetra::CombineMode::REPLACE);
      auto A_vv     = A.getMultiVectorView();
      updateJacobian(y, *Jac_, t);
      Jac_->apply(B_sm_vv, A_vv);
    }


    void applyJacobianFD(const state_t &U, const dense_matrix_type & A,scalar_t t, dense_matrix_type &JA)const{
      double eps = 1.e-5;
      state_t Up(U,Teuchos::DataAccess::Copy);
      velocity_type V0(*hyperMap_,3);
      velocity_type V_perturb(*hyperMap_,3);
      V_perturb.putScalar(0.);
      std::cout << "using FD jacobian" << std::endl;
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

    //auto dum = Tpetra_BlockMultiVector::makePointMap(*hyperMap_,3);
    auto hyperPointMapWithStencil = Tpetra::BlockMultiVector<scalar_t,lo_t,go_t>::makePointMap(*hyperMapWithStencil_,3);
    auto hyperPointMap = Tpetra::BlockMultiVector<scalar_t,lo_t,go_t>::makePointMap(*hyperMap_,3);
    using map_t = typename Tpetra::BlockMultiVector<scalar_t,lo_t,go_t>::map_type;
    const auto hyperPointMap_rcp = Teuchos::rcp(new map_t(hyperPointMap));
    const auto hyperPointMapWithStencil_rcp = Teuchos::rcp(new map_t(hyperPointMapWithStencil));
    //dense_matrix_type dum1(*hyperMap_,3,30);
    //dense_matrix_type dum2(*hyperMapWithStencil_,3,30);
    //auto pointMap1 = dum1.makePointMap(*hyperMap_,3);
    //auto pointMap2 = dum2.makePointMap(*hyperMapWithStencil_,3);
    //importer_ = std::make_shared<importer_t>( hyperMapWithStencil_, hyperMap_);
    importer_ = std::make_shared<importer_t>(hyperPointMapWithStencil_rcp, hyperPointMap_rcp);

    U_ = std::make_shared<nativeVec>(*hyperMapWithStencil_,3);

    xGrid_ = std::make_shared<Tpetra::Vector<>>(hyperMapWithStencil_);
    yGrid_ = std::make_shared<Tpetra::Vector<>>(hyperMapWithStencil_);
    for (int sid = 0; sid < sampleMeshPlusStencilSize_; sid++){
        auto gid = hyperMapWithStencil_->getGlobalElement(sid);
        auto ij = get_ij_from_gid(gid);
        xGrid_->replaceGlobalValue(gid, dx_*ij[0] + dx_*0.5);
        yGrid_->replaceGlobalValue(gid, dy_*ij[1] + dy_*0.5);
    }


    // Jacobian
    // construct a graph for the block matrix
    int nonZrPerRow_ = 5;
    crs_graph_type dataGraph(hyperMap_,nonZrPerRow_);
    assembleGraph(dataGraph);
    const lo_t blockSize = 3;
    Jac_ = std::make_shared<jacobian_type>(dataGraph,3);
    //Jac_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);


    };


  bool in_sample_mesh(go_t gid_in){
    bool val = false;
    const auto gids = hyperMap_->getMyGlobalIndices();
    for (int sid = 0; sid < sampleMeshSize_; sid++){
      auto gid = gids[sid];
      if (gid == gid_in){
        val = true;
        break;
      }
    }
    return val;
  }


  bool in_sample_mesh_plus_stencil(go_t gid_in){
    bool val = false;
    const auto gids = hyperMapWithStencil_->getMyGlobalIndices();
    for (int sid = 0; sid < sampleMeshPlusStencilSize_; sid++){
      auto gid = gids[sid];
      if (gid == gid_in){
        val = true;
        break;
      }
    }
    return val;
  }




  void assembleGraph(crs_graph_type & graph)
  {
   
    using tarr_it = Teuchos::ArrayView<go_t>;
    const auto gids = hyperMap_->getMyGlobalIndices();
    for (int sid=0; sid < sampleMeshSize_; sid++){
      std::vector<go_t> row_gids(0);
      auto gid = gids[sid];
      auto ij = get_ij_from_gid(gid);
      auto gid_im1 = get_gid_from_ij(ij[0] - 1,ij[1]);
      auto gid_ip1 = get_gid_from_ij(ij[0] + 1,ij[1]);
      auto gid_jm1 = get_gid_from_ij(ij[0],ij[1] - 1);
      auto gid_jp1 = get_gid_from_ij(ij[0],ij[1] + 1);
      row_gids.push_back(gid);
      go_t gid_sz = 1;
      if (in_sample_mesh(gid_im1)){
        row_gids.push_back(gid_im1);
        gid_sz += 1;
      }
      if (in_sample_mesh(gid_ip1)){
        row_gids.push_back(gid_ip1);
        gid_sz += 1;
      }
      if (in_sample_mesh(gid_jm1)){
        row_gids.push_back(gid_jm1);
        gid_sz += 1;
      }
      if (in_sample_mesh(gid_jp1)){
        row_gids.push_back(gid_jp1);
        gid_sz += 1;
      }
      
      //row_gids[0] = gid;
      //row_gids[1] = gid_im1;
      //row_gids[3] = gid_ip1;
      //row_gids[2] = gid_jm1;
      //row_gids[4] = gid_jp1;
      tarr_it tempArray(row_gids.data(),gid_sz);
      graph.insertGlobalIndices(gid, tempArray);
      //graph.insertLocalIndices(sid, tarr_it_lot(row_lids.data(),1));

    }
    
    graph.fillComplete();
  }//end







public:

  nativeVec getGaussianIC(double mu){
    auto xGridv = xGrid_->getData();
    auto yGridv = yGrid_->getData();
    int i = 0;
    for (int i=0; i < sampleMeshPlusStencilSize_; i++){
        double * uVal;
        U_->getLocalRowView(i,uVal);
        uVal[0] = 1. + mu*exp( - ( pow( xGridv[i] - 1.5,2)  + pow(yGridv[i] - 1.5,2) ));
        uVal[1] = 0.;
        uVal[2] = 0.;
    }
    return *U_;
  }

};

#endif
