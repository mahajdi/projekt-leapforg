#include<dune/pdelab/localoperator/idefault.hh>

#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/pattern.hh>

#include<dune/pdelab/finiteelement/localbasiscache.hh>


/** Lokalni operator za zadaću:
 *          un = uoul
 *          unn=uoldold
 *  u^{n+1} * phi + (dt)^2 * grad u^{n+1} * grad phi - (dt)^2 * f * phi - 2 * u^n * phi + u^{n-1} * phi = 0
 *
 *                  u = g   na \Gamma_D\partial\Omega
 *
 *
 *
 * \tparam BCType - klasa koja određuje tip rubnog uvjeta.
 */



template <typename BCType, typename FEM, class DGF>
class StationaryLocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume<StationaryLocalOperator<BCType,FEM, DGF> >,
  public Dune::PDELab::NumericalJacobianVolume<StationaryLocalOperator<BCType,FEM, DGF > >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<StationaryLocalOperator<BCType,FEM, DGF> >,
  public Dune::PDELab::NumericalJacobianBoundary<StationaryLocalOperator<BCType,FEM, DGF> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:

  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };

  StationaryLocalOperator(BCType& bctype_, // boundary cond.type
                         unsigned int intorder_=2, DGF& uold_, DGF& uoldold_ ) :
    bctype( bctype_ ), intorder( intorder_ ), uold (uold_), uoldold (uoldold_)
  {}

  // volume integral depending on test and ansatz functions
  // eg   = element (geometry) (klasa Dune::PDELab::ElementGeometry)
  // lfsu = lokalni prostor funkcija za rješenje (klasa Dune::PDELab::LocalFunctionSpace)
  // lfsv = lokalni prostor funkcija za test funkciju  (klasa Dune::PDELab::LocalFunctionSpace)
  // x    = vektor koeficijenata rješenja  (klasa Dune::PDELab::LocalVector)
  // r    = lokalni rezidual (određen "pogled" na  Dune::PDELab::LocalVector)
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    //  lfsu == lfsv
    const int dim = EG::Geometry::coorddimension;
    using Gradient = Dune::FieldVector<double,dim>;

    auto gt = eg.geometry().type();
    const auto & rule = Dune::QuadratureRules<double,dim>::rule(gt,intorder);

    for (const auto & ip : rule)
      {
        // Izračunaj bazne funkcije
        auto& phi = cache.evaluateFunction(ip.position(), lfsu.finiteElement().localBasis());

        // rješenje
        double u=0.0;   //sadašnji korak
        for (size_t i=0; i<lfsu.size(); ++i){
            u += x(lfsu,i)*phi[i];
        }


        // Gradijenti baznih funkcija
        auto& gradphihat = cache.evaluateJacobian(ip.position(), lfsu.finiteElement().localBasis());

        // grad g_K^{-t}
        const auto & jac = eg.geometry().jacobianInverseTransposed(ip.position());
        // Gradijenti baznih funkcija na fizičkom elementu.
        std::vector<Gradient> gradphi(lfsu.size());
        for (size_t i=0; i<lfsu.size(); i++)
          jac.mv(gradphihat[i][0],gradphi[i]);

        // Gradijent rješenja u
        Gradient gradu(0.0);
        for (size_t i=0; i<lfsu.size(); ++i)
          gradu.axpy(x(lfsu,i),gradphi[i]);

        // evaluate parameters;
        // auto globalpos = eg.geometry().global(ip.position());
        double f = 0;

        // integrate grad u * grad phi_i + a*u*phi_i - f phi_i
        double factor = ip.weight()*eg.geometry().integrationElement(ip.position());
        double un;   //prošli korak
        double unn;  //pretprošlikorak
        uold.evaluate(eg,ip.position(), un );
        uoldold.evaluate(eg,ip.position(), unn );
        for (size_t i=0; i<lfsv.size(); ++i)
          r.accumulate(lfsv, i, (u*phi[i] + 0.01*gradu*gradphi[i] - 0.01 * f*phi[i]-2*un*phi[i] + unn*phi[i] ) * factor);
      }
  }
/*
  void preStep (double time, double dt, int stages) {
    bctype.setTime(time+dt); // ako treba promijeni Dirichletovu granicu
    Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::preStep(time,dt,stages);
  }*/
private:
  BCType& bctype;
  unsigned int intorder;
  typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
  DGF& uold;
  DGF& uoldold;
};
