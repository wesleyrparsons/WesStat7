unit CLI;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Bivariate,
  Correlations,
  DataDisplay,
  DataManager,
  Globals,
  LinearAlgebra,
  Math,
  Multivariate,
  Nominal,
  Parser,
  PCA,
  Regression,
  SysUtils,
  WesUnit,
  Univariate;

procedure RunCLI;

implementation

var
  SlashMark: Integer;                                 // Counts the Tokens; pinpoints the /.
  Line, Token: String;                                // The Line, Token, and error message.
  i, Pos: Integer;                                    // Pos traverses the Line.
  AddError, Okay, IsVar, HitOptions: Boolean;
  InArguments: Boolean;                               // Flag true if in arguments before SlashMark + 1, the objects;
  EvalLine, Col: WVectorType;
  EigenValues: CVectorType;
  WQ1, WQ2, WQ3, WL1, WL2, WL3, WL4: FloatType;
  Lower, Upper, AInv, Q, RR, C: WMatrixType;
  U, S, VT, EigenVectors: CMatrixType;
  OPrecision, OWidth, OMLIter: Integer;
  OSample, OExt, OUseIntercept: Boolean;
  OMissingData: MissingDataType;
  OMatInvEpsilon, ORegConvEpsilon, ORankTolerance: FloatType;
  r, adjr, cov, v1, v2: FloatType;
  Cond: FloatType;
  rLower90, rUpper90, rLower95, rUpper95, rLower99, rUpper99: FloatType;
  // Error messages.
  ErrMess: Array[0..14] of String = (' is not an option', ' is not a command',                    //  0  1
    ' is not a variable', ': One variable and two numbers are required.',                         //  2  3
    ': Variable or trim values are not correct.', 'One variable and one number are required.',    //  4  5
    'File name is required.', 'Two variables are required.',                                      //  6  7
    'Extra arguments are ignored.', ': is not a constant.',                                       //  8  9
    ': Two variables are required.', ': One variable and a % value are required.',                // 10 11
    ': One variable and an integer k value are required.', ': Row is not found.',                 // 12 13
    'Input is required.');                                                                        // 14

function FoundVar(const x: TokenArgType): Boolean;
begin
  Result := False;
  if x.V >= 0 then Result := True;
end;

// Tells if a S is a scalar name, and returns the index number.
function FoundScalIndex(const x: TokenArgType; var Index: Integer): Boolean;
var
  i: Integer;
begin
  Result := False;
  for i:= 0 to nScal - 1 do
    if UpCase(SData[i].Name) = UpCase(x.S) then begin
      Result := True;
      Index := i;
      Break
    end;
end;

function FoundScal(const x: TokenArgType): Boolean;
begin
  if x.I <= nScal then
    Result := True
  else
    Result := False;
end;

function FoundNumber(const x: TokenArgType): Boolean;
begin
  if x.I > 0 then
    Result := True
  else
    Result := False;
end;

//For ease of use, put WToken variables into WVar and Col.
procedure PutTokensIntoWVar;
var
  i: Integer;
Begin
  //if not FoundVar(WTokenArg[1]) then Exit;
  if nArg < 1 then Exit;
  SetLength(WVar, nArg);
  for i := 1 to nArg do                // nArg = 1 is first non-command token (argument).
    WVar[i - 1] := WTokenArg[i].V;     // WVar is zero based.
  nVar := nArg;
end;

// For options, which can function as options or as commands.
procedure OptionsCase;
begin
  HitOptions := True;
  Case UpCase(WTokenArg[iOpt].S) of
    // For options that cannot take arguments. If Main Case (iOpt = 0) and WTokenArgs exist, then show error.
    'NULL': writeln('Null does nothing.');
    'DISPLAYWVAR', 'WVAR' : DisplayWVar;
    'SIGNIFICANCE', 'SIG': begin ShowSignificance := True; writeln('Significance: ', ShowSignificance); if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]); end;
    'NOSIGNIFICANCE', 'NOSIG': begin ShowSignificance := False; writeln('Significance: ', ShowSignificance); if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]); end;
    'INTERCEPT', 'INT': begin UseIntercept := True; writeln('Use intercept: ', UseIntercept); if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]); end;
    'NOINTERCEPT', 'NOINT': begin UseIntercept := False; writeln('Use intercept: ', UseIntercept); if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]); end;
    'EXTENDED', 'EXT': begin Ext := True; writeln('Statistics extended.'); end;
    'NOEXTENDED', 'NOEXT': begin Ext := False; writeln('Statistics not extended.'); end;
    'SAMPLE', 'SAM': begin Sample := True; writeln('Sample: ', Sample); if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]); end;
    'POPULATION', 'POP': begin Sample := False; writeln('Sample: ', Sample); if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]); end;
    'MDMEAN': begin MissingData := MDMean; ConvertNaNinData; if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]); end;  // Add writeln to show state.
    'MDMEDIAN', 'MDMED': begin MissingData := MDMedian; ConvertNaNinData; if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]); end;
    'MDZERO', 'MD0': begin MissingData := MDZero; ConvertNaNinData; if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]); end;

    // For options that take one number argument. If Main Case (iOpt = 0) and WTokenArgs > 1, then show error.
    'PRECISION', 'PREC': if WTokenArg[iOpt + 1].I > 0 then begin
        Precision := WTokenArg[iOpt + 1].I;
        writeln('Precision = ', Precision);
        {if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln('Precision = ', Precision);     If no argument, write out prior setting.}
      end
      else writeln('Positive integer required (default is 5).');
   'WIDTH': if WTokenArg[iOpt + 1].I > 0 then begin
        Width := WTokenArg[iOpt + 1].I;
        writeln('Width = ', Width);
        //if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]);
      end
      else writeln('Positive integer required (default is 14).');
    'RANKTOLERANCE', 'RANKTOL': if WTokenArg[iOpt + 1].F > 0 then begin
        RankTolerance := WTokenArg[iOpt + 1].F;
        writeln('Rank tolerance: ', RankTolerance);
        //if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]);
      end
      else writeln('Small positive number required (default is 1e-9).');
    'MATRIXINVERSEEPSILON', 'MATINVEPS', 'INVEPS':  if WTokenArg[iOpt + 1].F > 0 then begin
        MatInvEpsilon := WTokenArg[iOpt + 1].F;
        writeln('Matrix inversion epsilon: ', MatInvEpsilon);
        //if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]);
      end
      else writeln('Small positive number required (default is 1e-10).');
    'REGRESSIONCONVERGENCEEPSILON', 'REGEPS': if WTokenArg[iOpt + 1].F > 0 then begin
        RegConvEpsilon := WTokenArg[iOpt + 1].F;
        writeln('Regression convergence epsilon: ', RegConvEpsilon);
        //if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]);
      end
      else writeln('Small positive number required (default is 1e-9).');
    'REGRESSIONITERATIONS', 'REGITER': if WTokenArg[iOpt + 1].I > 0 then begin
        MLIter := WTokenArg[iOpt + 1].I;
        //if not (WTokenArg[1].S = EmptyStr) and (iOpt = 0) then writeln(ErrMess[8]);
        writeln('Regression iterations: ', MLIter);
      end
      else writeln('Positive integer required(default is 500).');
    else HitOptions := False;
  end;
  Inc(iOpt);  // So loop iterates to next WToken in options.
end;

// For commands.
procedure MainCase;
var
  i: Integer;
begin
  // For linear algebra commands that cannot take arguments.
  Case UpCase(WTokenArg[0].S) of
    'INVERSE', 'INV':  begin
       InvertWMatrixType(WData, AInv, Okay);
       if Okay then
         DisplayMatrixM('Inverse ', Ainv, T2);
       if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
    end;
    'LUDECOMPOSITION', 'LUDECOMP', 'LU': begin
      LUDecomposition(WData, Lower, Upper);
      writeln('Matrix L');
      DisplayMatrixM('', Lower, T2);
      writeln('Matrix U');
      DisplayMatrixM('', Upper, T2); pause;
      if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
    end;
    'QRDECOMPOSITION', 'QRDECOMP','QR': begin
      QRDecomposition(WData, Q, RR);
      writeln('Matrix Q');
      DisplayMatrixM('', Q, T2);
      writeln('Matrix R');
      DisplayMatrixM('', RR, T2);
      if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
    end;
    'CHOLESKY', 'CHO': begin
      CholeskyDecomposition(WData, C);
      if C <> nil then begin
        writeln('Matrix L-Transpose');
        DisplayMatrixM('', C, T2);
      end;
      if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
    end;
    'SINGULARVALUEDECOMPOSITION', 'SVD': begin
      SVD(WData, U, S, VT);
      writeln('Matrix Sigma');
      DisplayMatrixM('', U, T2);
      writeln('Matrix S');
      DisplayMatrixM('', S, T2);
      writeln('Matrix V-Transpose');
      DisplayMatrixM('', VT, T2);
      if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
    end;
    'EIGENDECOMPOSITION', 'JACOBI','EIGEN': begin
      JacobiEigenDecomposition(CreateATA(WData), EigenVectors, EigenValues);
      ReportEigenResults(EigenVectors, EigenValues);
      if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
    end;
    'SPACES', 'SP': begin Spaces(WData); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'PRINCIPALCOMPONENT', 'PRINCIPALCOMPONENTCORR', 'PCA', 'PCACORR': begin
      PCAnalysis(CorrPCA);
      if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
    end;
    'PRINCIPALCOMPONENTCOV', 'PCACOV': begin
      PCAnalysis(CovPCA);
      if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
    end;
    'DETERMINANT', 'DET': begin writeln('The determinant is ', Determinant(WData), '.'); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'CONDITION', 'COND': begin
      ConditionNumber(WData, Cond, Okay);
      if not Okay then writeln('Condition number not computable.')
      else writeln('The condition number is ', SmartFloat(Cond));
      if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
    end;
    'FROBENIUS', 'FRO': begin writeln('The Frobenius Norm is ', FrobeniusNormW(WData), '.'); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'MATRIX', 'MAT': begin MatrixFacts(WData); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'TRACE', 'TR': begin writeln('The trace is ', Trace(WData), '.'); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'RANK': begin writeln('The rank is ', MatrixRank(WData), '.'); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'SYMMETRIC', 'SYM': begin
      write('The matrix is ');
      if IsSymmetric(WData) = False then write('not ');
      writeln('symmetric.');
      if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
     end;
    'SQUARE', 'SQ': begin
      write('The matrix is ');
      if IsSquare(WData) = False then write('not ');
      writeln('square.');
      if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]);
     end;

    // For other commands that cannot take arguments.
    'MANUALINPUT','MANINPUT': ManuallyInputData;
    'FILLCONST', 'FILLCONSTANT': begin FillScalars; if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'PCONST', 'PCONSTANT': begin DisplayAllScalars(T2); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'PVAR', 'PVARIABLE' : for i := 0 to nCol - 1 do DisplayVectorM(WData[i].Name, WData[i], T2);
    'ERROR', 'ERR': begin DisplayErrorStackM; if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'I', 'INF': begin ReportProgramState; if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'D': begin DebugOn := not DebugOn; writeln('Debug is '); if DebugOn then writeln('on.') else write ('off.'); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'V': begin Verbose := not Verbose; writeln('Verbose is '); if Verbose then writeln('on.') else write ('off.'); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'VERBOSE': begin Verbose := True; if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'NOVERBOSE', 'NOV': begin Verbose := False; if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'OPTIONS', 'OPT', 'OPTS', 'SETTING': begin ShowOptions; if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'P1': begin DisplayMatrixM('Main Data', WData, T1); DisplayAllScalars(T1);  if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'P', 'P2': begin DisplayMatrixM('Main Data', WData, T2); DisplayAllScalars(T2);  if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'P3': begin DisplayMatrixM('Main Data', WData, T3); DisplayAllScalars(T3);  if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'C1': begin DisplayAllScalars(T1); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'C', 'C2': begin DisplayAllScalars(T2); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'C3': begin DisplayAllScalars(T3); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'MAINDATA', 'MAIN': begin DisplayMatrixM('Main Data ', WData, T2); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'T1' : begin GetTestData(WTokenArg[0].S); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'T2' : begin GetTestData(WTokenArg[0].S); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'T3' : begin GetTestData(WTokenArg[0].S); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'T4' : begin GetTestData(WTokenArg[0].S); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'T5' : begin GetTestData(WTokenArg[0].S); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'T6' : begin GetTestData(WTokenArg[0].S); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'T7' : begin GetTestData(WTokenArg[0].S); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'T8' : begin GetTestData(WTokenArg[0].S); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'T9' : begin GetTestData(WTokenArg[0].S); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;
    'T10' : begin GetTestData(WTokenArg[0].S); if not (WTokenArg[1].S = EmptyStr) then writeln(ErrMess[8]); end;

    // For regression commands, which take more than one argument.
    'OLSREGRESSION', 'REGRESSION', 'OLSREG', 'OLS': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
        RegressionMode := MLEOLS;                                              // Add type of regression as parameter.
        MainRegression;
      end
    else writeln(WTokenArg[i].S, ErrMess[7]);
    'MLEREGRESSIONGD','GRADIENTDESCENT',
    'MLEGD', 'GD': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
        RegressionMode := MLEGD;
        MainRegression;
    end
    else writeln(WTokenArg[i].S, ErrMess[7]);
    'MLEREGRESSIONNR', 'NEWTONRAPHSON',
    'MLENR', 'NR': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
        RegressionMode := MLENR;
        MainRegression;
      end
      else writeln(WTokenArg[i].S, ErrMess[7]);
    'MLEREGRESSIONNEST', 'NESTEROV',
    'MLENEST', 'NEST': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
        RegressionMode := MLENest;
        MainRegression;
      end else writeln(WTokenArg[i].S, ErrMess[7]);

    // For functions that take one variable (and here can be multiple variables in sequence).
    'NOMINAL', 'NOM': if FoundVar(WTokenArg[i]) then NominalStatistics else writeln(WTokenArg[i].S, ErrMess[1]);  // Using WVar.
    'UNIVARIATE', 'U', 'UNI': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then DisplayUnivariateStatistics(WTokenArg[i]) else writeln(WTokenArg[i].S, ErrMess[1]);

    // The following use WVar and nVar.
    'CORRELATIONMATRIX', 'PEARSONMATRIX', 'RMAT', 'CORRELATION', 'CORR': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      MultivariateCorrelations(Pearson);
    end
    else writeln(WTokenArg[i].S, ErrMess[1]);
    'SPEARMANMATRIX', 'SMAT': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      MultivariateCorrelations(Spearman);
    end
    else writeln(WTokenArg[i].S, ErrMess[1]);
    'KENDALLMATRIX', 'KMAT': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      MultivariateCorrelations(Kendall);
    end
    else writeln(WTokenArg[i].S, ErrMess[1]);
    'HOEFFDINGMATRIX', 'HMAT': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      MultivariateCorrelations(Hoeffding);
    end
    else writeln(WTokenArg[i].S, ErrMess[1]);
    'MUTUALINFORMATIONMATRIX', 'MIMAT': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      MultivariateCorrelations(MI);
    end
    else writeln(WTokenArg[i].S, ErrMess[1]);
    'NORMALIZEDMUTUALINFORMATIONMATRIX', 'NMIMAT': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      MultivariateCorrelations(NMI);
    end
    else writeln(WTokenArg[i].S, ErrMess[1]);
    'COVARIANCEMATRIX', 'COVMAT': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      MultivariateCorrelations(Covar);
    end
    else writeln(WTokenArg[i].S, ErrMess[1]);

    { Measures of Central Tendency }
    'MULTIVARIATE', 'MULTI', 'MUL', 'M': begin
      for i := 1 to nArg do if not FoundVar(WTokenArg[i]) then begin writeln(WTokenArg[i].S, ErrMess[2]); Exit; end;
      MultivariateStatistics;
    end;
    'MEAN', 'ARITHMETICMEAN', 'AMEAN': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Arithmetic Mean = ', SmartFloat(ArithmeticMean(WData[WTokenArg[i].V])), '  ')
      else writeln(WTokenArg[i].S, ErrMess[2]);
    'GEOMETRICMEAN', 'GMEAN': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].I].Name, ' Geometric Mean = ', SmartFloat(GeometricMean(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'HARMONICMEAN', 'HMEAN': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Harmonic Mean = ', SmartFloat(HarmonicMean(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'MEANPAIRWISEDISTANCE', 'MPAIR': for i := 1 to nArg do  if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Mean Pairwise Distance = ', SmartFloat(MeanPairwiseDistance(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'ROOTMEANSQUARE', 'RMS': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Root Mean Square = ', SmartFloat(RootMeanSquare(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'CUBICMEAN', 'CUBIC': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Cubic Mean = ', SmartFloat(CubicMean(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'COEFFICIENTVARIATION', 'COEFFVAR': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then writeln(WData[WTokenArg[i].V].Name, ' Coefficient Of Variation = ', SmartFloat(CoefficientOfVariation(WData[WTokenArg[i].V], Sample)), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'MEANABSOLAUTEDEVIATION', 'MAD': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Mean Absolaute Deviation = ', SmartFloat(MeanAbsoluteDeviation(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'STANDARDERRORMEAN', 'STDERRMEAN': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Standard Error of Mean = ', SmartFloat(StandardErrorOfMean(WData[WTokenArg[i].V], Sample)), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);

    {Measures of Dispersion}
    'VARIANCE': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Variance = ', SmartFloat(WesUnit.Variance(WData[WTokenArg[i].V], Sample)), '  ');
    'STANDARDDEVIATION', 'STDDEV', 'SD': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Std. Dev. = ', SmartFloat(WesUnit.StandardDeviation(WData[WTokenArg[i].V], Sample)), '  ');
    'MEDCOUPLE', 'MEDC': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Medcouple = ', SmartFloat(Medcouple(WData[WTokenArg[i].V])), '  ');
    'SKEWNESS', 'SKEW': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Skewness = ', SmartFloat(Skewness(WData[WTokenArg[i].V], Sample)), '  ');
    'KURTOSIS', 'KURT': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Kurtosis = ', SmartFloat(Kurtosis(WData[WTokenArg[i].V], Sample)), '  ');
    'EXCESSKURTOSIS', 'EXKURT': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Excess Kurtosis = ', SmartFloat(ExcessKurtosis(WData[WTokenArg[i].V], Sample)), '  ');
    'BIMODALITYCOEFFICIENT', 'BIMOD': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Bimodality Coefficient = ', SmartFloat(BimodalityCoefficient(WData[WTokenArg[i].V])), '  ');
    'SYMMETRYINDEX', 'SYMIND', 'SI': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Symmetry Index = ', SmartFloat(SymmetryIndex(WData[WTokenArg[i].V])), '  ');
    'HILLESTIMATOR', 'HILLEST', 'HE': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Hill Estimator = ', SmartFloat(HillEstimator(WData[WTokenArg[i].V], Trunc(Sqrt(nRow)))), '  ');

    {Measures of Equality}
    'GINICOEFFICIENT', 'GINI' : for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Gini Coef. = ', SmartFloat(GiniCoefficient(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'LORENZMUENZNERCOEFFICIENT', 'LORMUNZ' : for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Lorenz Muenzner Coef. = ', SmartFloat(LorenzMuenzner(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'THEILINDEX', 'THEIL' : for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Theil Index = ', SmartFloat(TheilIndex(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);

    {NON-PARAMETRIC STATISTICS}
    {Measures of Central Tendency}
    'MEDIAN', 'MED': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Median = ', SmartFloat(Median(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'MEDIANPAIRWISEDISTANCE', 'MEDPAIR': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Median Pairwise Distance = ', SmartFloat(MedianPairwiseDistance(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'HIGHMEDIAN', 'HIMED': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' High Median = ', SmartFloat(HighMedian(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'LOWHMEDIAN', 'LOMED': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Low Median = ', SmartFloat(LowMedian(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'TRIMEDIAN', 'TRIMED': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Trimedian = ', SmartFloat(Trimedian(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'MODE': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Mode = ', SmartFloat(Mode(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'HUBERESTIMATOR', 'HUBER': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Huber Estimator = ', SmartFloat(HuberEstimator(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'HODGESLEHMANNESTIMATOR', 'HL': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Hodges Lehmann Estimator = ', SmartFloat(HodgesLehmann(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'TUKEYBIWEIGHT', 'TUKEYBW': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Tukey Biweight = ', SmartFloat(TukeyBiweight(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);

    {Measures of Range}
    'MINIMUM', 'MIN': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Minimum = ', SmartFloat(MinmValue(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'MAXIMUM', 'MAX': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Maximum = ', SmartFloat(MaxmValue(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'RANGE', 'RAN': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Range = ', SmartFloat(WesUnit.Range(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'MIDRANGE', 'MIDRAN': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Midrange = ', SmartFloat(MidRange(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'FIRSTQUARTILE', 'FIRSTQ': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
       DeepCopyVector(WData[WTokenArg[i].V], Col);
       Quartiles(Col, WQ1, WQ2, WQ3);
       writeln(WData[WTokenArg[i].V].Name, ' First Quartile = ', SmartFloat(WQ1), '  ');
    end
    else writeln(WTokenArg[i].S, ErrMess[2]);
    'THIRDQUARTILE', 'THIRDQ': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
       DeepCopyVector(WData[WTokenArg[i].V], Col);
       Quartiles(Col, WQ1, WQ2, WQ3);
       writeln(WData[WTokenArg[i].V].Name, ' Third Quartile = ', SmartFloat(WQ3), '  ');
     end
     else writeln(WTokenArg[i].S, ErrMess[2]);
    'MIDHINGE', 'MIDH': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
       DeepCopyVector(WData[WTokenArg[i].V], Col);
       Quartiles(Col, WQ1, WQ2, WQ3) ;
       writeln(WData[WTokenArg[i].V].Name, ' Midhinge = ', SmartFloat(Midhinge(WQ1, WQ3)), '  ');
     end
     else writeln(WTokenArg[i].S, ErrMess[2]);
    'Q3/Q1RATIO', 'Q31R': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
       DeepCopyVector(WData[WTokenArg[i].V], Col);
       Quartiles(Col, WQ1, WQ2, WQ3);
       writeln(WData[WTokenArg[i].V].Name, ' Q3/Q1 Ratio = ', SmartFloat(Q3Q1Ratio(WQ1, WQ3)), '  ');
    end
    else writeln(WTokenArg[i].S, ErrMess[2]);
    'INTERQUARTILERANGE', 'INTERQR' : for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' InterQuartile Range = ', SmartFloat(InterquartileRange(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'SEMIINTERQUARTILERANGE', 'SEMIINTERQR': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
       DeepCopyVector(WData[WTokenArg[i].V], Col);
       Quartiles(Col, WQ1, WQ2, WQ3);
       writeln(WData[WTokenArg[i].V].Name, ' Semi-Interquartile Range = ', SmartFloat(SemiInterquartileRange(WQ1, WQ3)), '  ');
     end
     else writeln(WTokenArg[i].S, ErrMess[2]);
    'P90-P10RATIO', 'P9010RATIO': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' P90-P10 Ratio = ', SmartFloat(P90P10Ratio(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'P95-P5RATIO', 'P955RATIO': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' P95-P5 Ratio = ', SmartFloat(P95P5Ratio(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'P75-P25SPREAD', 'P7525S': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' P75-P25 Spread = ', SmartFloat(P75P25Spread(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'P95-P5SPREAD', 'P955S': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' P95-P5 Spread = ', SmartFloat(P95P5Spread(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'RANGERATIO', 'RANGER': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Range Ratio = ', SmartFloat(RangeRatio(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'QUARTILECOFFICIENTOFDISPERSION', 'QCOEFFDISP': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
       DeepCopyVector(WData[WTokenArg[i].V], Col);
       Quartiles(Col, WQ1, WQ2, WQ3);
       writeln(WData[WTokenArg[i].V].Name, ' Quartile Coef. of Dispersion = ', SmartFloat(QuartileCoefficientOfDispersion(WQ1, WQ3)), '  ');
    end
    else writeln(WTokenArg[i].S, ErrMess[2]);
    'YULESKEWNESS', 'YULESKEW': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      DeepCopyVector(WData[WTokenArg[i].V], Col);
      Quartiles(Col, WQ1, WQ2, WQ3);
      writeln(WData[WTokenArg[i].V].Name, ' Yule Skewness = ', SmartFloat(YuleSkewness(WData[WTokenArg[i].V], WQ1, WQ3)), '  ');
    end
    else writeln(WTokenArg[i].S, ErrMess[2]);
    'SNESTIMATOR', 'SNEST': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Rousseeuw Croux Sn Estimator = ', SmartFloat(SnEstimator(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'QNESTIMATOR', 'QNEST': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Rousseeuw Croux Qn Estimator = ', SmartFloat(QnEstimator(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'MEANABSOLUTEDEVIATIONFROMMEAN', 'MABSDEVM': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Mean Absolute Deviation From Mean = ', SmartFloat(MeanAbsoluteDeviationFromMean(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'MEDIANABSOLUTEDEVIATIONFROMMEAN', 'MEDABSDEVM': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Median Absolute Deviation From Mean = ', SmartFloat(MedianAbsoluteDeviationFromMean(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'MEANABSOLUTEDEVIATIONFROMMEDIAN', 'MABSDEVMED': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Mean Absolute Deviation From Median = ', SmartFloat(MeanAbsoluteDeviationFromMedian(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'MEDIANABSOLUTEDEVIATIONFROMMEDIAN', 'MEDABSDEVMED': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Median Absolute Deviation From Median = ', SmartFloat(MedianAbsoluteDeviationFromMedian(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'STANDARDERROROFMEDIAN', 'STDERRMED': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Standard Error of Median = ', SmartFloat(StandardErrorOfMedian(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'SPEARMANRANKSKEWNESS', 'SRS': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Spearman''s Rank Skewness = ', SmartFloat(SpearmanRankSkewness(WData[WTokenArg[i].V])), '  ') else writeln(WTokenArg[i].S, ErrMess[2]);
    'BOWLEYSKEWNESS', 'BOWLEY': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      DeepCopyVector(WData[WTokenArg[i].V], Col);
      Quartiles(Col, WQ1, WQ2, WQ3);
      writeln(WData[WTokenArg[i].V].Name, ' Bowley''s Skewness = ', SmartFloat(BowleySkewness(WQ1, WQ2, WQ3)), '  ');
    end;
    'KELLYSKEWNESS', 'KELLY': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      writeln(WData[WTokenArg[i].V].Name, ' Kelly''s Skewness = ', SmartFloat(KellySkewness(WData[WTokenArg[i].V])), '  ');
    'LSCALE': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      DeepCopyVector(WData[WTokenArg[i].V], Col);
      Quartiles(Col, WQ1, WQ2, WQ3);
      LMoments(Col, WL1, WL2, WL3, WL4);
      writeln(WData[WTokenArg[i].V].Name, ' L-Scale = ', SmartFloat(WL2), '  ');
     end
     else writeln(WTokenArg[i].S, ErrMess[2]);
    'LSKEWNESS', 'LSKEW': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      DeepCopyVector(WData[WTokenArg[i].V], Col);
      Quartiles(Col, WQ1, WQ2, WQ3);
      LMoments(Col, WL1, WL2, WL3, WL4);
      writeln(WData[WTokenArg[i].V].Name, ' L-Skewness = ', SmartFloat(WL3), '  ');
     end
     else writeln(WTokenArg[i].S, ErrMess[2]);
    'LKURTOSIS', 'LKURT': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      DeepCopyVector(WData[WTokenArg[i].V], Col);
      Quartiles(Col, WQ1, WQ2, WQ3);
      LMoments(Col, WL1, WL2, WL3, WL4);
      writeln(WData[WTokenArg[i].V].Name, ' L-Kurtosis = ', SmartFloat(WL4), '  ');
     end
     else writeln(WTokenArg[i].S, ErrMess[2]);
    'LMOMENTS': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      DeepCopyVector(WData[WTokenArg[i].V], Col);
      Quartiles(Col, WQ1, WQ2, WQ3);
      LMoments(Col, WL1, WL2, WL3, WL4);
      writeln(WData[WTokenArg[i].V].Name, ' L-Scale = ', SmartFloat(WL2), '   L-Skewness = ',
        SmartFloat(WL3), '   L-Kurtosis = ', SmartFloat(WL4), '  ');
    end
     else writeln(WTokenArg[i].S, ErrMess[2]);

    { Utility Functions }
    'CENTER', 'CENT': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      CenterColumn(WData[WTokenArg[i].V]);
      writeln('Data have been centered.');
     end
     else writeln(WTokenArg[i].S, ErrMess[2]);
    'STANDARDIZE', 'STAND': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then begin
      StandardizeColumn(WData[WTokenArg[i].V]);
      writeln('Data have been normalized.');
     end
     else writeln(WTokenArg[i].S, ErrMess[2]);
     'DELVARIABLE', 'DELVECTOR', 'DELCOLUMN', 'DELVAR', 'DELCOL', 'DELVEC': if FoundVar(WTokenArg[i]) then     // For dels, what if name of var is in S
     for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
        DeleteVariable(WTokenArg[i].V) else writeln(WTokenArg[i].S, ErrMess[2]);
    'DELCONST', 'DELCONSTANT': for i := 1 to nArg do if FoundScalIndex(WTokenArg[i], WTokenArg[i].V) then DeleteScalar(WTokenArg[i].V) else writeln(WTokenArg[i].S, ErrMess[10]);
    'DELROW', 'DELOBS': for i := 1 to nArg do if WTokenArg[i].I <= nRow then DeleteRow(WTokenArg[i].I - 1) else writeln(WTokenArg[i].I, ErrMess[13]);
    'COL1', 'VECTOR1', 'VEC1': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      DisplayVectorM(WData[WTokenArg[i].V].Name, WData[WTokenArg[i].V], T1) else writeln(WTokenArg[i].S, ErrMess[2]);
    'COL', 'VECTOR', 'VEC', 'COL2', 'VECTOR2', 'VEC2' : for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      DisplayVectorM(WData[WTokenArg[i].V].Name, WData[WTokenArg[i].V], T2) else writeln(WTokenArg[i].S, ErrMess[2]);
    'COL3', 'VECTOR3', 'VEC3': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then
      DisplayVectorM(WData[WTokenArg[i].V].Name, WData[WTokenArg[i].V], T3) else writeln(WTokenArg[i].S, ErrMess[2]);
    'CHECKCOLUMN', 'CHKCOL': for i := 1 to nArg do if FoundVar(WTokenArg[i]) then CheckColumn('Check Column', WData[WTokenArg[i].V]) else writeln(WTokenArg[i].S, ErrMess[2]);
    'Q', 'QUIT', 'X', 'EXIT': Exit;

    // For command with two variable arguments.
    'BIVARIATE', 'B', 'BIV': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then
      BivariateCorrelations(WTokenArg[1].V, WTokenArg[2].V)
      else writeln(WTokenArg[i].S, ErrMess[10]);


    'PEARSON', 'PEARSONSR', 'CORRCOEFF': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := PearsonsR(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      if (nRow > 2) then
        adjr := 1 - ((1 - r * r) * (nRow - 1) / (nRow - 2))
      else
        adjr := 0.0;
      writeln('Pearson''s R = ', r:2 :Precision);
      writeln('   Significance =  ', PearsonsRSignificance(r, Length(WData[WTokenArg[1].V].Value)):2 :Precision);
      writeln('Pearson''s Adjusted R = ', adjr:2 :Precision);
      writeln('   Significance =  ', PearsonsRSignificance(adjr, Length(WData[WTokenArg[1].V].Value)):2 :Precision);
      writeln('Pearson''s R-Squared = ', (r * r):2 :Precision);
      writeln('   Pearson''s Adjusted R-Squared = ', (adjr * adjr):2 :Precision);
      // Compute confidence intervals
      PearsonCI(r, nRow, 0.90, rLower90, rUpper90);
      PearsonCI(r, nRow, 0.95, rLower95, rUpper95);
      PearsonCI(r, nRow, 0.99, rLower99, rUpper99);
      // Display confidence intervals
      writeln('Confidence intervals at 90%.     Upper: ', SmartFloat(rUpper90), '   Lower: ', SmartFloat(rLower90));
      writeln('Confidence intervals at 95%.     Upper: ', SmartFloat(rUpper95), '   Lower: ', SmartFloat(rLower95));
      writeln('Confidence intervals at 99%.     Upper: ', SmartFloat(rUpper99), '   Lower: ', SmartFloat(rLower99));
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'COVARIANCE', 'COV': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      v1 := Variance(WData[WTokenArg[1].V], Sample);
      v2 := Variance(WData[WTokenArg[2].V], Sample);
      cov := Covariance(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      writeln('Variance of first variable: ', SmartFloat(v1));
      writeln('Variance of second variable: ', SmartFloat(v2));
      writeln('Covariance of variables: ', SmartFloat(cov));
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'SPEARMAN', 'SPEARMANSRHO', 'SPEAR': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := SpearmansRho(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      writeln('Spearman''s Rho = ', r:2 :Precision);
      writeln('   Significance =  ', SpearmansRhoSignificance(WData[WTokenArg[1].V], WData[WTokenArg[2].V], r):2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'KENDALL', 'KENDALLSTAU', 'KENDTAU': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := KendallsTau(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      writeln('Kendall''s Tau = ', r:2 :Precision);
      writeln('   Significance =  ', KendallsTauSignificance(WData[WTokenArg[1].V], WData[WTokenArg[2].V]):2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'HOEFFDINGSD', 'HOEFFDING', 'HOEFF': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := HoeffdingD(WData[WTokenArg[1].V], WData[WTokenArg[2].V], 1000, 500);
      writeln('Hoeffding''s D = ', r:2 :Precision);
      writeln('   Significance =  ', HoeffdingsDSignificance(WData[WTokenArg[1].V], WData[WTokenArg[2].V], 1000, 500, 1000, r):2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'MUTUALINFORMATION', 'MUTUAL', 'MI': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := MutualInformation(WData[WTokenArg[1].V], WData[WTokenArg[2].V], Trunc(Sqrt(Length(WData[WTokenArg[1].V].Value))));
      Writeln('Mutual Information = ', r:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'NORMALIZEDMUTUALINFORMATION', 'NORMMUTUAL', 'NMI': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := NormalizedMutualInformation(WData[WTokenArg[1].V], WData[WTokenArg[2].V], Trunc(Sqrt(Length(WData[WTokenArg[1].V].Value))));
      Writeln('Normalized Mutual Information = ', r:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'DISTANCECORRELATION', 'DISTCORR', 'DC': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := DistanceCorrelation(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      Writeln('Distance Correlation = ', r:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'JENSENSHANNONDIVERGENCE', 'JENDIV', 'JSDIV': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      v1 := JensenShannonDivergence(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      Writeln('Jensen Shannon Divergence = ', V1:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'JENSENSHANNONDISTANCE', 'JENDIST', 'JSDIST': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      v2 := JensenShannonDistance(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      Writeln('Jensen Shannon Distance = ', V2:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'KULLBACKLEIBERDIVERGENCE', 'KULLDIV', 'KLDIV': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := KLDivergence(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      Writeln('Kullback Leibler Divergence = ', r:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'HELLINGERDISTANCE', 'HELLDIST', 'HDIST': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      R := HellingerDistance(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      Writeln('Hellinger Distance = ', r:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'TOTALVARIATIONDISTANCE', 'TOVARIST', 'TVDIST': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := TotalVariationDistance(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      Writeln('Total Variation Distance = ', r:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'BHATTACHARYYACOEFFICIENT', 'BHATTCOEFF', 'BC': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := BhattacharyyaCoefficient(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      Writeln('Bhattarchayya Coefficient = ', r:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'BHATTACHARYYADISTANCE', 'BHATTDIST', 'BD': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := BhattacharyyaDistance(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      Writeln('Bhattarchayya Distance = ', r:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'CHISQUAREDIVERGENCE', 'CHIDIV', 'CHID': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := ChiSquareDivergence(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      Writeln('Chi Square Divergence = ', r:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'CANBERRADISTANCE', 'CANBDIST', 'CD': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      r := CanberraDistance(WData[WTokenArg[1].V], WData[WTokenArg[2].V]);
      Writeln('Canberra Distance = ', r:2 :Precision);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'MANNWHITNEYU', 'MWU': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      iVar1 := WVar[0];
      iVar2 := WVar[1];
      MannWhitneyU(WData[iVar1], WData[iVar2]);
    end else writeln(WTokenArg[i].S, ErrMess[10]);
    'SIGNTEST', 'ST': if FoundVar(WTokenArg[1]) and FoundVar(WTokenArg[2]) then begin
      iVar1 := WVar[0];
      iVar2 := WVar[1];
      SignTest(WData[iVar1], WData[iVar2]);
    end else writeln(WTokenArg[i].S, ErrMess[10]);

    // For commands with two number arguments.
    'CREATE', 'CR' : if (WTokenArg[1].I > 0) and (WTokenArg[2].I > 0) then begin
      nCol := WTokenArg[1].I;
      nRow := WTokenArg[2].I;
      CreateRandomData;
    end;
    'WINSORIZE', 'WIN': if FoundVar(WTokenArg[1]) and (WTokenArg[2].F > 0.0) and (WTokenArg[3].F > 0.0) then begin  // Check args.
         WinsorizeWColumn(WData[WTokenArg[1].V], WTokenArg[2].F, WTokenArg[3].F);
         writeln('Data have been been winsorized.');
       end
     else writeln(WTokenArg[0].S, ErrMess[3]);
    'TRIM': if (WTokenArg[1].I <= nRow - 1) and FoundVar(WTokenArg[2]) and FoundVar(WTokenArg[3]) then
       DataTrim(WTokenArg[1].I, WTokenArg[2].V, WTokenArg[3].V)
     else writeln(ErrMess[4]);

    // For commands with one string argument.
    'READ': if not (WTokenArg[1].S = EmptyStr) then ReadCSVData(WTokenArg[1].S, nRow, nCol) else writeln(ErrMess[6]);
    'SAVE': if not (WTokenArg[1].S = EmptyStr) then SaveData(WTokenArg[1].S) else writeln(ErrMess[6]);
    'PATH': if not (WTokenArg[1].S = EmptyStr) then begin FilePath := WTokenArg[1].S; writeln('Path = ', FilePath); end else writeln(ErrMess[14]);
    'LOG': if not (WTokenArg[1].S = EmptyStr) then LogtoFile(True, WTokenArg[1].S) else writeln(ErrMess[6]);
    'NOLOG': if not (WTokenArg[1].S = EmptyStr) then LogtoFile(False, WTokenArg[1].S) else writeln(ErrMess[6]);
    'MISSINGDATA', 'MD': Case WTokenArg[1].S of
      'MEAN': begin MissingData := MDMean; ConvertNaNinData; end;
      'MEDIAN': begin MissingData := MDMedian; ConvertNaNinData; end;
      'ZERO', '0': begin MissingData := MDZero; ConvertNaNinData; end;
      else writeln('Missing data option not found.');
    end;

    // For commands with one variable argument, and one number argument.
    'YEOJOHNSON', 'YJ', 'YEO': if FoundVar(WTokenArg[1]) and (WTokenArg[2].I > 0) then begin
        YeoJohnsonTransform(WData[1], WTokenArg[2].F);                         // Ck args.
        writeln('Data have been been transformed by Yeo-Johnson.');
    end else writeln(ErrMess[5]);
    'SHANNONENTROPY', 'SHAN' : if FoundVar(WTokenArg[1]) and FoundNumber(WTokenArg[2]) then
      writeln(WData[WTokenArg[i].V].Name, ' Shannon Entropy = ', SmartFloat(ShannonEntropy(WData[WTokenArg[1].V], WTokenArg[2].I)), '  ')
    else writeln(ErrMess[12]);
    'GENERALIZEDMEAN', 'GENMEAN': if FoundVar(WTokenArg[1]) and FoundNumber(WTokenArg[2]) then
      writeln(WData[WTokenArg[0].V].Name, ' Generalized Mean with d = ', WTokenArg[2].F: 4: 2, ': ',
        SmartFloat(GeneralizedMean(WData[WVar[0]], WArg[nArg - 1])))
    else writeln(ErrMess[6]);
    'QUANTILE', 'QUANT': if FoundVar(WTokenArg[1]) and FoundNumber(WTokenArg[2]) then
      writeln(WData[WTokenArg[0].V].Name, ' Quantile with p = ', SmartFloat(WVar[nArg - 1]), ': ',
        SmartFloat(Quantile(WData[WTokenArg[1].V], WTokenArg[2].V)))
    else writeln(ErrMess[6]);
    'PERCENTILE', 'PERCENT': for i := 1 to (WTokenArg[i].I - 1) div 2 do begin          // fix
      if FoundVar(WTokenArg[i]) and (WTokenArg[i].I < (2 * i + 2)) then begin
        writeln(WData[WTokenArg[i].V].Name, ' Percentile with % = ', WTokenArg[i + 1].I, ': ',
          SmartFloat(Percentile(WData[WTokenArg[i].V], WTokenArg[i + 1].F / 100)), '  ');
      end else if not HitOptions then writeln(ErrMess[6]);
    end;
    'HARRISDAVISQUANTILE', 'HDQ': for i := 1 to (WTokenArg[i].I - 1) div 2 do begin      // fix
      if FoundVar(WTokenArg[i]) and (WTokenArg[i].I < (2 * i + 2)) then begin
        writeln(WData[WTokenArg[i].V].Name, ' Harris-Davis Quantile with p = ', WTokenArg[i + 1].I, ': ',
          SmartFloat(HarrisDavisQuantile(WData[WTokenArg[i].V], WTokenArg[i + 1].F)), '  ');
      end else writeln(ErrMess[6]);
    end;
    'PEARSONKMOMENT', 'PKM':  if FoundVar(WTokenArg[1]) and (WTokenArg[2].I > 0) then
      writeln(WData[WTokenArg[0].V].Name, ' Pearson''s Kth Moment with k = ', Trunc(WTokenArg[2].F), ' is ',
        SmartFloat(PearsonKthMoment(WData[WVar[0]], Trunc(WTokenArg[2].F))))
    else writeln(ErrMess[12]);
    'QUICKSELECT', 'QS': if FoundVar(WTokenArg[1]) and (WTokenArg[2].I > 0) then
        writeln(WData[WTokenArg[0].V].Name, ' Kth Largest Element = ', SmartFloat(WTokenArg[0].F))
      else writeln(WTokenArg[i].S, ErrMess[12]);

    // For no command found.
    else if not HitOptions then
      writeln('Command not found.');
  end;
end;

function IsAlphaNum(const C: Char): Boolean;
begin
  Result := IsLetter(C) or IsDigit(C) or (C in ['-', '+', '.']);
end;

procedure ArchiveOptions;
begin
  OPrecision := Precision;
  OWidth := Width;
  OSample := Sample;
  OExt := Ext;
  OMissingData := MissingData;
  OUseIntercept := UseIntercept;
  ORankTolerance := RankTolerance;
  OMatInvEpsilon := MatInvEpsilon;
  ORegConvEpsilon := RegConvEpsilon;
  OMLIter := MLIter;
end;

procedure RestoreOptions;
begin
  Precision := OPrecision;
  Width := OWidth;
  Sample := OSample;
  Ext := OExt;
  MissingData := OMissingData;
  UseIntercept := OUseIntercept;
  RankTolerance := ORankTolerance;
  MatInvEpsilon := OMatInvEpsilon;
  RegConvEpsilon := ORegConvEpsilon;
  MLiter := OMLiter;
end;

// Report the state of WTokenArg.
procedure ReportAllTokenState;
var
  i: Integer;
begin
  writeln('WTokenArg Full Report: iTA = ', iTA, '; nArg = ', nArg, '; Line Position = ', Pos, '; Length of Line = ', Length(Line), '; SlashMark = ', SlashMark,
    '; iArg = ', iArg, '; iOpt = ', iOpt, '; nOpt = ', nOpt);
  i := iTA - 1;
  for i := 0 to iTA - 1 do
    writeln('iTA=', i, '; nArg=', nArg - 1,   '; SlashMark=', SlashMark,
    '.  String=', WTokenArg[i].S, '; Integer=', WTokenArg[i].I,
      '; Float=', WTokenArg[i].F: 4: 2, '; Variable=', WTokenArg[i].V);
  DisplayWVar;
end;

procedure ReportCurrentTokenState;
var
  i: Integer;
begin
  i := iTA;
  writeln('Token="', Token, '"  iTA=', i, '; nArg=', nArg - 1,   '; SlashMark=', SlashMark,
  '.  String=', WTokenArg[i].S, '; Integer=', WTokenArg[i].I,
    '; Float=', WTokenArg[i].F: 4: 2, '; Variable=', WTokenArg[i].V);
end;

// Get the next Token.
procedure GetTokenArg;
begin
  // Step thru non-Alphanum characters.
  while (Pos <> Length(Line) + 1) and not (IsAlphaNum(Line[Pos])) do begin
    // If signal /, then mark it to show end of arguments and start of options.
    if Line[Pos] = '/' then begin
      If InArguments then SlashMark := iTA - 1;   // Zero based.
      InArguments := False;                       // If not in arguments part, then second slash.
    end;
    Inc(Pos);
  end;

  // Get alphanumeric characters and build Token.
  while (Pos < Length(Line) + 1) and (Line[Pos] <> '/') and (Line[Pos] <> ' ')
    and (Line[Pos] <> '=') do begin    // Include = because of options.
    // Get a character from Line and add to Token.
    if IsAlphaNum(Line[Pos]) then
      Token := Token + Line[Pos];
    Inc(Pos);
    // At this point, Pos is one more than the last char of Token.
  end;
  if DebugOn then
    ReportCurrentTokenState;
end;

// Put the Token into WTokenArg.
procedure PutTokenArg;
var
  i, n, Code: Integer;
begin
  // Initialize WTokenArg.
  WTokenArg[iTA].F := NaN;                       // Start with F not being a number.
  WTokenArg[iTA].I := -99;                       // Same with I, 99 is a sentinel.
  WTokenArg[iTA].S := Token;                     // Put Token into String.
  if iTA = 0 then Exit;                          // If Command, just get the command into S.
  Val(Token, n, Code);                           // Get n, a float if Token is a number.
  WTokenArg[iTA].I := Trunc(n);                  // Then put n into Integer.  Code doesn't matter.
  WTokenArg[iTA].V := -99;                       // -99 is sentinel to show unassigned.
  WTokenArg[iTA].F := n;                         // Then put n into Float.

  // First, see if any token is an integer in the nCol range, and InArguments.
  for i := 0 to nCol - 1 do
    if (n > 0) and (n <= nCol) then begin        // Is it in range?  N is one based, so subtract 1.
      WTokenArg[iTA].V := n - 1;                 // Put it in V.
      WTokenArg[iTA].S := WData[n - 1].Name;     // Put its name in S.
      Exit;                                      // Don't put it in I or V, may be variable, like create 12 33.
    end;

  // If not integer in range, then, second, see if any token has a name in the nCol range, and inArgument.
  if InArguments then
    for i := 0 to nCol - 1 do
      if UpCase(WData[i].Name) = UpCase(Token) then begin
        WTokenArg[iTA].V := i;
        Exit;
      end;
  end;

// If Token = 'ALL'.
procedure PutAllTokenArg;
var
  i: Integer;
begin
  for i := 0 to nCol - 1 do begin
    Token := WData[i].Name;
    PutTokenArg;
    Inc(iTA);
    Inc(nArg);
  end;
end;

// Declare a variable if Command is DCL etc.
function DeclareVariable: Boolean;
var
  Command: String;
begin
//  Command := UpCase(WTokenArg[iTA].S);
  Command := UpCase(Token);
  if (Command = 'DCL') or (Command = 'DECLARE') or (Command = 'DCLV') or (Command = 'VAR')
    or (Command = 'VARIABLE') then begin
      Result := True;
      IsVar := True;                   // Keep track if Token is a variable or a constant.
    end
  else begin
    Result := False;
    IsVar := False;
  end;
end;

// Declare a scalar if Command is Const etc.
function DeclareScalar: Boolean;
var
  Command: String;
begin
  Command := UpCase(WTokenArg[iTA].S);
  if (Command = 'CONST') or (Command = 'CON') or (Command = 'CONSTANT') or (Command = 'DCLC')
    or (Command = 'C') or (Command = 'DCLC') then begin
      Result := True;
      IsVar := False;                  // Keep track if Token is a variable or a constant.
    end
  else begin
    Result := False;
    IsVar := False;
  end;
end;
// 0           1        2              3        4         5          6            7         8            9
// mean        a        size           3        ht /      sig        prec         11        nosig        sig
// WTA[0]=mean WTA[1]=a WTA[2]=size    WTA[3]=3 WTA[4]=ht WTA[5]=sig WTA[6]=prec  WTA[7]=11 WTA[8]=nosig WTA[9]=sig
// Command = mean       nArg = 3 (0-3) SlashMark = 5      iAT = 10   iOpt = 6     nOpt = 10
// Args go from 1 to 4; Opts go from 6 to 10
// Main program for this unit.
procedure RunCLI;
begin
  ReportDataState;           // Reports whether data exists.

  while true do begin        // Loop to read each Line.  Outer loop.
    ArchiveOptions;          // Save the old options.
    for i := 0 to MaxToken do
      WTokenArg[i].S := EmptyStr;
    iTA := 0;
    InArguments := True;
    SlashMark := 0;
    nArg := -1;
    nOpt := 0;
    writeln;
    Write('> ');             // no no two loops one to read line one to process tokens

    // Read a Line.
    ReadLn(Line);
    Line := Trim(Line);      // Trim off spaces.

    // Easy exit out of Line loop for Tokens.
    if (Line = 'Q') or (Line = 'QUIT') or (Line = 'X') or (Line = 'EXIT') then begin
      //Facade := MDI;
      Exit;
    end;

    // Easy restart if Line is not right.
    if (Line = EmptyStr) or (Line = '/') or not IsAlphaNum(Line[1]) then Continue;

    // Pos is the parse position in the Line. Initialize Pos and SlashMark.
    Pos := 1;

    repeat                   // Loop to read each Token and put in WTokenArg. Inner Loop.

      Token := EmptyStr;
      // Get a TokenArg;
      GetTokenArg;

      // If Command is Var, Dcl, Const etc. then take a detour to Parserr.
      if DeclareVariable or DeclareScalar then begin
        if WData = nil then begin                          // Necessary because if variable added, must know number of rows.
          write('No data exists. Please enter number of rows: ');
          readln(nRow);
        end;

        Delete(Line, 1, Pos);                              // Remove first chars that are DCL etc.
        //Line := Copy(Line, 1, Pos('/', Line + '/') - 1); // Remove chars after any slash. Do I need this?
        if DebugOn then
          writeln('Sending to ParserExpression: *', Line, '*');
        EvalLine := ParserExpression(Line);                 // Get result vector from Parserr.
        if DebugOn then begin
          writeln('Returning from ParserExpression. EvalLine: ', Length(EvalLine.Value), ' value ');
          DisplayVectorM(EvalLine.Name, EvalLine, T3);
        end;                                               // End detour to Parserr, but need to add variable or scalar.

        // Check Name is good and add variable or scalar.
        if (EvalLine.Name <> '') and IsProperVariableName(EvalLine.Name)
          and IsProperScalarName(EvalLine.Name)then begin
          if IsVar then
            AddWVector(EvalLine.Name, EvalLine, AddError)
          else
            AddWConstant(EvalLine.Name, EvalLine, AddError);

          // State the variable or constant has or has not been added.
          if not AddError then begin
            write(EvalLine.Name, ' declared as ');
            If IsVar then
              writeln('variable.')
            else
              writeln('constant.');
            if Verbose or DebugOn then begin
              writeln('After adding variable or constant.');
              if IsVar then
                DisplayVectorM(EvalLine.Name, EvalLine, T3)
              else
                AddWConstant(EvalLine.Name, EvalLine, AddError);
            end;
            Continue;         // Skip rest of while loop. Or continue to Options? Test whether this works.
          end
          else begin
            writeln('Variable or constant not declared.');
            Continue;        // Skip loop.
          end;
        end;
      end;                   // End adding variable or scalar from Parserr.

      // Get Token, and add to WTokenArg.
      //Inc(Pos);            // Not sure Inc is right.
      if UpCase(Token) = 'ALL' then begin
        PutAllTokenArg;      // Put all the variable into WTokenArg.
      //Inc(iTA, nCol);
      end
      else begin
        PutTokenArg;         // Put one argument there.
        Inc(iTA);
        if InArguments then Inc(nArg) else Inc(nOpt);
      end;

      // Report WTokenArg -- everything done above.
    until Pos >= Length(Line);         // End the loop thru the Tokens.

    // Execute options in cases -- do before executing Command.
    if (SlashMark > 0) and (SlashMark < iTA) then begin
      iOpt := SlashMark + 1;           // Start with the first option at iTA + 1.
      while iOpt <= iTA - 1 do          // Loop until end of TokenArg.
        OptionsCase;
    end;

    // Preprocess WData for NaN data.
    ConvertNaNinData;

    // Execute the main cases.
{    if SlashMark = 0 then    // If SlashMark not found, that is, zero, then put SlashMark at end of iTA.
      SlashMark := iTA - 1;}
    HitOptions := True;      // Reversed if else statement shows hit.
    //nTA := SlashMark;      // nTA is number of WTokenArgs, same as SlashMark, less one.
    //iTA := 0;              // iTA loops thru the WTokenArgs. Wait, does it? i loops.
    //nArg := iTA - 1;
    iOpt := 0;               // iOpt loops thru the Options.
    nOpt := iTA - SlashMark;
    PutTokensIntoWVar;       // Helpful for routines that use WVar, Col, etc.
    if DebugOn then begin
      ReportAllTokenState;
      Pause;
    end;

    // Run Options Case and then Main Case.
    OptionsCase;
    MainCase;
    if SlashMark > 0 then    // Not working well.
      RestoreOptions;        // Restore the original options. But only if options exist.}

  end;                       // End the while loop thru the Line.
end;

end.

