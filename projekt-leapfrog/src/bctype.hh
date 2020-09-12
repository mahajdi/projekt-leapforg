#include <dune/pdelab/constraints/conforming.hh>

/*  Selekcija Dirichletove granice */
class BCTypeParam
  : public Dune::PDELab::DirichletConstraintsParameters
{
    // Dirichletova granica se može micati u vremenu.
  double time;
public:
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
                   ) const
  {
    //auto xg = intersection.geometry().global( coord );
/*
    if( xg[0]>1.0-1E-6 )
      return false; // x=1 nije Dirichletova granica */
    return true;     // sve je Dirichletova granica
  }

  //! Postavi vrijeme.
  void setTime (double t) { time = t; }

};
