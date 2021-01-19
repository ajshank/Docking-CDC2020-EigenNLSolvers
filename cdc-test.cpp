/*
  Eigen/C++ implementation for solving the non-linear
  algebraic equations presented in:
    "Dynamic Path Generation for Multirotor Aerial Docking in Forward Flight"
    Conference on Decision & Control (CDC), 2020.
    
  There are six non-linear equations with six unknowns; together, these
  form a parametric representation for optimal docking paths.
  
  Compile as:
  $> g++ cdc-test.cpp -o cdctest $(pkg-config --cflags eigen3) -O3 -ffast-math

  -- aj // Nimbus Lab, 2020.
*/

#include <iostream>

#include <chrono>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>

using namespace Eigen;
typedef std::chrono::high_resolution_clock ClockTime;
typedef std::chrono::time_point<ClockTime> ClockTimePoint;
typedef std::chrono::microseconds uSeconds;

// functor base-class
template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
struct TFunctor
{
  typedef _Scalar Scalar;
  enum
  {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  const int m_inputs, m_values;

  TFunctor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime)
  {}
  TFunctor(int inputs, int values) : m_inputs(inputs), m_values(values)
  {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

};

struct cdc_functor : TFunctor<double, 6, 6>
{
    double c1, c2, c3, c4, t1, t2;
    double sq2;
    double F, L0, R, w1, w2;
    void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
    
    cdc_functor(void) : TFunctor<double, 6, 6>()
    {
      sq2 = std::sqrt(2);
      // problem constants -- can load through constructor
      F = 1.0;
      L0 = 10.0;
      R = 2.0;
      w1 = 3.0;
      w2 = 5.0;
    }

    int operator()(const VectorXd &x, VectorXd &fvec)
    {
        c1 = x(0); c2 = x(1); c3 = x(2); c4 = x(3);
        t1 = x(4); t2 = x(5);
        // eqn 9--14
        fvec <<   c1+c2 - (F-L0+R),
                  c1*exp(-t1/sq2) + c2*exp(t1/sq2),
                  c3*exp(-t1/sq2) + c4*exp(t1/sq2),
                  c3*exp(-t2/sq2) + c4*exp(t2/sq2) - R,
                  (c2*c2 - c4*c4)*exp(t1*sq2) - 0.5*(w1-w2),
                  2*c4*exp(t2/sq2) - ( R + sqrt(2*w2) );
        
        return 0;
    }
    int df(const VectorXd &x, MatrixXd &J)
    {
      c1 = x(0); c2 = x(1); c3 = x(2); c4 = x(3);
      t1 = x(4); t2 = x(5);
       
      J <<
        1, 1, 0, 0, 0, 0,
	      exp(-t1/sq2), exp(t1/sq2), 0, 0, -(sq2*exp(-t1/sq2)*(c1 - c2*exp(sq2*t1)))/2,  0,
	      0, 0, exp(-t1/sq2), exp(t1/sq2), -(sq2*exp(-t1/sq2)*(c3 - c4*exp(sq2*t1)))/2,  0,
	      0, 0, exp(-t2/sq2), exp(t2/sq2), 0, -(sq2*exp(-t2/sq2)*(c3 - c4*exp(sq2*t2)))/2,
	      0, 2*c2*exp(sq2*t1), 0, -2*c4*exp(sq2*t1), sq2*exp(sq2*t1)*(c2*c2 - c4*c4),  0,
        0,  0,   0, 2*exp(t2/sq2), 0,  sq2*c4*exp((sq2*t2)/2);
        
      return 0;
    }
};

int main()
{
  const int n=6;
  int info;
  VectorXd x(n);

  // initial guess
  x << -1, 1, -1, 1, 2, 2;

  // dogleg: trust region
  cdc_functor prob1;
  
  HybridNonLinearSolver<cdc_functor> tr(prob1);
  tr.parameters.factor = 1.0;
  tr.diag.setConstant(n, 1.0);
  tr.useExternalScaling = true;
    
  ClockTimePoint ts, te;
  ts = ClockTime::now();
    info = tr.solve(x);
  te = ClockTime::now();
  
  std::cout.precision(4); std::cout.setf(std::ios::fixed);
  double dt = std::chrono::duration_cast<uSeconds>(te - ts).count();
  std::cout << "solver status: " << info << std::endl;
  std::cout << "solution:\t" << x.transpose() << std::endl;
  std::cout << "Minimized F:\t" << tr.fvec.transpose() << std::endl;
  std::cout << "f-eval, j-eval:\t" << tr.nfev << ", " << tr.njev << std::endl;
  std::cout << "solver (us) " << dt << std::endl << "-- --\n";
  
  // LM
  // initial guess
  x << -1, 1, -1, 1, 2, 2;
  LevenbergMarquardt<cdc_functor> lm(prob1);
  info = lm.minimize(x);

  std::cout << "solver status: " << info << std::endl;
  std::cout << "solution:\t" << x.transpose() << std::endl;
  std::cout << "Minimized F:\t" << lm.fvec.transpose() << std::endl;
  std::cout << "f-eval, j-eval:\t" << lm.nfev << ", " << lm.njev << std::endl;

  return 0;
  
}
