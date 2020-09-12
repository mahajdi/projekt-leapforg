#include <cmath>
#include <dune/pdelab/common/function.hh>


// egzaktno rješenje
// uexact(x,y,t)=x^2-y^2+t
template <int dim>
double exact(Dune::FieldVector<double, dim> const & glob, double time_t){

    double rez = 0.0;
    const double x = glob[0];
    const double y = glob[1];
    rez = x*x -y*y+ time_t;
    return rez;

}



// Klasa obilježja. Ovo je samo pokrata koja čini kod čitljivijim.
template<typename GV>
using ScalarTraits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>;

/*  Dirichletov rubni uvjet proširen na čitavu domenu.
 *  U t=0 daje inicijalni uvjet!
 */
template<typename GV>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<ScalarTraits<GV>, BCExtension<GV> >
{
  const GV& gv;
  double  time;
public:
  using Traits = ScalarTraits<GV>;

  // Sačuvaj gridview
  BCExtension (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    auto  x= e.geometry().global(xlocal);
   /* if (x[0]<1E-6 && x[1]>0.25-1e-6 && x[1]<0.5+1e-6)
      y = std::sin(time);
    else
      y = 0.0; */
    //if (time == 0)
    y=exact(x,time);
    //else y= 0.0;
    return;
  }

  // Referenca na gridview
  inline const GV& getGridView () {return gv;}

  // Postavljanje vremena. Potrebno pozvati prije evaluate() !
  void setTime (double t) {time = t;}
};


/*
template <typename GV>
using ATraits = Dune::PDELab::AnalyticGridFunctionTraits<GV,double,1>;// 1 = skalarna funkcija

// Klasa koja daje diskretnu mrežnu funkciju egzaktnog rješenja.
template <typename GV>
class ExactGF : public Dune::PDELab::AnalyticGridFunctionBase<ATraits<GV>, ExactGF<GV>>
{
    double time;
public:
   typedef Dune::PDELab::AnalyticGridFunctionBase<ATraits<GV>, ExactGF<GV> > BaseT;

   ExactGF(GV const & gv) : BaseT(gv) {}

   void evaluateGlobal(typename ATraits<GV>::DomainType const & x, typename ATraits<GV>::RangeType & y) const
   {
       y = exact(x,time);
   }
   // Bazna klasa daje metodu evaluate().
    void setTime(double t) {time =t;}
};

*/

// za početnu egzaktnu brzinu
template <typename GV>
class ExactGFv : public Dune::PDELab::AnalyticGridFunctionBase
             <
               Dune::PDELab::AnalyticGridFunctionTraits<GV,double,1>, // 1 = skalarna funkcija
               ExactGFv<GV>
             >
{
public:
   typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,double,1> Traits;
   typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, ExactGFv<GV> > BaseT;

   ExactGFv(GV const & gv) : BaseT(gv) {}
   void evaluateGlobal(typename Traits::DomainType const & x, typename Traits::RangeType & y) const
   {
       //y = DiffusionParameter<GV>::exact(x);
       y = 0.1*1;
   }
};

