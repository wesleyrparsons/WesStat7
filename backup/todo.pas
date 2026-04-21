unit ToDo;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  Crt;

implementation

begin

  { Changes High Priority }
  // Finish Parser recode.
  // In Options section, special settings for iterations and tolerance, and check epsilons and iterations throughout.
  // In CLI, pass bt tabs, #9.     while (Pos <= Length(Line)) and (Line[Pos] in [' ', #9]) do Inc(Pos);
  // Change procs to functions for vectors and matrices. Change names to nouns and verbs.

  { Changes High Priority }
  // Linear Algebra:
  //   Express inverse matrix as determinant times inverse; add adjugate. T10, e.g.
  //   Add spectral norm (largest singular value))
  //   Positive definite / semidefinite check; Orthogonal / unitary check; Idempotent check (A² = A); Nilpotent check (Aᵏ = 0 for some k).
  //   For inverting matrices, do I use functions or procedures?
  //   I have a matrix inversion routine just for checking. Fix.
  // Univariate & Bivariate statistics.
  //   Add reps to Hoeffding's D.
  //   Correlations: check and update, very carefully, the significance functions.
  //   Test univariate functions with small datasets, and 1 obs and 1 var.
  //   Check WPMoment.
  //   Add QS to advanced uni stats. Same for Hill estimator. Same for Entropy with bins.
  //   I have two entropy routines. One can take a second parameter and could be advanced.
  //   Put clean up functions in the statistics function itself, rather than in Univariate and Bivariate.
  //   Outlier procs need to be checked more.
  //   Bimodality Coefficient. Dip/Silverman are formal tests: stronger evidence. Entropy/divergence measures give shape context.
  //     Clustering metrics reveal latent structure. BC pairs beautifully with: Dip Test, Silverman’s Test, Hellinger / JSD,
  //     Shannon entropy, and medcouple-based skewness diagnostics.
  // Regression.
  //   Add logistic and probit regression.
  //   Update routines: Momentum, Nesterov Acceleration, Adam (simplest modern optimizer), Coordinate Descent,
  //     Barzilai–Borwein Step, BFGS, and Levenberg–Marquardt.
  //   Add Options toggle for regression, as listwise deletion or mean/median replacement. Need code for missing data -- NaN?
  // Command Line and Parser.
  //   In CLI, rewrite using CLVar, CLInt, CLFloat and CLOpt. New repeat loop for clauses with ';'.
  //   Add dot as function in parser.
  // Others.
  //   Use AI to create a help file for uni stats.
  //   Add KNN toolkit, and use KNN for missing data.
  //   In read data, add check for vars in first row, and option to force on or off.
  //   Add Verbose and DebugOn options in LA and Regression and elsewhere where needed.
  //   Set eps tol and iter. // Add low and high eps from Correlations?

  { Guidelines }
  // Use ClrScr on tables: Univariate, Bivariate, Multivariate, Correlations.
  // Should I display a message in the interpreter to indicate an option has been set, such as prec or scinot? Yes.
  // Set a path in user options and interpreter to change directory path. It's FilePath: String.
  // I have an entropy function in Correlations.
  // Should I use IsEffectivelyEqual more often? Yes.
  // Use guards on all routines with matrices. Don't use infinity checks.
  // Use Nan cgecks only on significance functions, such as t, gamma, beta, CDF, erf.
  // Use AddToErrorStack and Result := NaN or nil. Remember Exit(Nan) and Exit(nil).
  // Use nouns for funcs and verbs for procs.
  // Why can't I get more than 25 lines? Use console options on console.

  { Will not do }
  // On Univariates, use fewer pauses. Nope. One pause is good.
  // Perhaps I should make every matrix a WMatrix. Maybe a function to generate a name, and a flag to show whether stats included?
  // No. The classic variables are WData -- other results are CMatriux and CVector.
  // The conversions happen when a WVector or WMatrix is multiplied by another matrix, like XtX or XtY.
  // Check Medcouple is operating correctly.  Yes, redone 11/28/2025.
  // Is format 8: 10 right in correlations? Yes, looks good. Do not allow overwrite by DecPrec in tables.
  // Maybe have a sort routine on simple vectors. Nope, not needed so far.
  // Don't do Mahalanobis or robust distances; output is vector.
  // Create vat NameWVar so name of var can be passed for corr cov matrices. Already done.
  // Change WesUnit and Univariate and Bivariate procs to use CMatrixType results. Nope.
  // Consider using exceptions. Nope.
  // Don't do weighted statistics; use evaluation instead.


 {                       1st Review    2nd Review     Wikipedia      Buxbaum      LLM      Eps/Iter/Nan/CU
  Univariate                 X                                          X          X               X
  Univariate Advanced        X                                          X          X               X
  Bivariate                  X                                                     X               X
  Multivariate               X                                                     X               X
  Nominal                    X                                                     X               X
  Correlations
  Precision
  Display Data
  Linear Algebra
  Regression
  PCA
  BasicFuncs
  OLS
  PCA
  Small Data (0/1)
  Weird Data
  Modify Data
 }
{Constant         Single         Double         Extended        Notes

Epsilon           1e-6           1e-12          1e-18
MachineEpsilon    1.192e-7       2.225e-16      1.084e-19
CleanUpEpsilon    1e-6           1e-12          1e-18           Equals Epsilon
Tolerance         1e-6           1e-9           1e-16
VeryTinyNumber    1e-40          1e-300         1e-4900
MinRealNumber     1.4016e-45     2.225e-308     3.645e-4951

Proposed Const
  MachineEpsilon            1.192e-7       2.225e-16      1.084e-19
  CompilerEpsilon           1e-6           1e-12          1e-18
  CompilerHighTolerance     1e-6           1e-12          1e-18
  CompilerTolerance         1e-6           1e-9           1e-16
  CompilerLowTolerance      1e-4           1e-6           1e-12

Proposed Var
  Epsilon           1e-6           1e-12          1e-18
  HighTolerance     1e-6           1e-12          1e-18
  Tolerance         1e-6           1e-9           1e-16
  LowTolerance      1e-4           1e-6           1e-12
  ?RankTolerance                    1e-6
  ?MatInvTolerance                  1e-12
  UserTolerance

  FastMaxIter       50;    // Newton, refinement
  MediumMaxIter     300;   // LM, QR, eigen
  SlowMaxIter       2000;  // gradient / EM
  UserMaxIter
Function              Uses
Effectively           Epsilon
ShannonEntropy        Epsilon15
MedianOfKernel        1e-12
SymmetryIndex         1e-12
HuberEstimator        Epsilon9         Iter1000
ChiSquare             Epsilon
Corr                  LowEpsilon = Epsilon6; HighEpsilon = Epsilon12;
IncompleteGamma       LowEpsilon       Iter100
IncBetaCF             HighEpsilon      Iter200
KendallsTau           1e-10
KendallsTauSign       1e-10
ComputeMidRanbks      HighEpsilon
ComputeJointRanbks    HighEpsilon
GetTieSums            HighEpsilon
CreateCorrFromCov     HighEpsilon
RandomPermuteTest                      Iter100
Reg
  ConditionHessianR   Epsilon6, 9
  CheckConvergence                                       Tol 1e-10
  MLERegression                        Iter5000
UpdateGD              1e-12
LA
  IsSymmetric                                            Tolerance
  RankFromR                                              RankTolerance
  NullSpace                                              RankTolerance
  QRDecomposition     Epsilon
  Cholesky            MatInvEpsilon
  Jacobi              MatInvEpsilon    Iter100
  Determinant         MatInvEpsilon
  MatrixRank          MatInvEpsilon, 1e-10
  Invert_Matrix       Epsilon15, 1e-12
  TestInvert          1e-12

MaxIter10000      10000
MaxIter1000       1000
MaxIter500        500
MaxIter200        200
MaxIter100        100
}

{ WesStat Proposed Data Types
                        Boolean	        Integer	        Float          String
  Simple Array Type			                Add C?
    None	        BType/Boolean	IType/Integer	FType/Float    String/Stype;
    Vector		                IVectorType	FVectorType    SVectorType;
    Array		                IMatrixType	FMatrixType

  Wes Record Array Type
    None         	BType/Boolean	IType/Integer	FType/Float
    Vector		BVectorType     WIVectorType	WFVectorType
    Array		                WIMatrixType	WFMatrixType
}



