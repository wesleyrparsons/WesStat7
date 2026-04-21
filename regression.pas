unit Regression;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Correlations,
  DataDisplay,
  DataManager,
  Globals,
  LinearAlgebra,
  Math,
  SysUtils,
  WesUnit;

procedure MainRegression(const WVar: IVectorType);

implementation

procedure SaveRegressionResults(var Residuals, YHat: WVectorType);
var
  AddError: Boolean;
begin
// Save residuals.
  write('Save residuals as additional variable? (y/n) ');
  readln(WInput);
  AddError := True;
  if UpCase(WInput) = 'Y' then begin
    write('Input Name: ');
    readln(WInput);
    Residuals.Name := WInput;
    // Check proper name -- or put this in Add proc
    AddWVector(WInput, Residuals, AddError);
    if not AddError then
      writeln('Variable ', WInput, ' has been added.')
    else
     writeln('Error: Variable ', WInput, ' was not added.')
  end;

  // Save YHat.
  write('Save y-hat as additional variable? (y/n) ');
  readln(WInput);
  AddError := True;
  if UpCase(WInput) = 'Y' then begin
    write('Input Name: ');
    readln(WInput);
    // Check proper name
    YHat.Name := WInput;
    AddWVector(WInput, YHat, AddError);
    if not AddError then
      writeln('Variable ', WInput, ' has been added.')
    else
     writeln('Error: Variable ', WInput, ' was not added.')
  end;
end;

function CheckConvergence(const OldBetas, NewBetas: CVectorType;
  OldLogSigma, NewLogSigma: FloatType; Tol: FloatType): Boolean;
var
  i, k: Integer;
  Diff, Rel, MaxDiff: FloatType;
begin
  k := Length(NewBetas);

  // Safety: empty vectors → treat as converged.
  if k = 0 then begin
    Result := True;
    Exit;
  end;

  MaxDiff := 0.0;
  // β convergence.
  for i := 0 to k - 1 do begin
    Diff := Abs(NewBetas[i] - OldBetas[i]);
    Rel  := Diff / (1.0 + Abs(OldBetas[i]));
    if Rel > MaxDiff then
      Maxdiff := Rel;
  end;

  // σ convergence via log σ.
  Diff := Abs(NewLogSigma - OldLogSigma);
  Rel  := Diff / (1.0 + Abs(OldLogSigma));
  if Rel > MaxDiff then
    Maxdiff := Rel;
  Result := (Maxdiff < Tol);
end;

function LogLikelihood(const Y: WVectorType; const Yhat: WVectorType; const Sigma: FloatType): FloatType;
var
  j, n: Integer;
  Resid, SS: FloatType;
begin
  n := Length(Y.Value);

  // Sum of squares.
  SS := 0.0;
  for j := 0 to n - 1 do begin
    Resid := Y.Value[j] - Yhat.Value[j];
    SS := SS + Resid * Resid;
  end;

  // Gaussian log-likelihood (constant term included).
  if Sigma <= 0 then
    Result := VeryTinyNumber  // Extremely low likelihood if invalid.
  else
    Result := -0.5 * n * Ln(TwoPi) - n * Ln(Sigma) - SS / (2.0 * Sqr(Sigma));
end;

procedure ComputePrediction(const X: WMatrixType; const Betas: CVectorType; out Yhat: WVectorType);
var
  i, j, n, k: Integer;
  sum: FloatType;
begin
  k := Length(X);
  n := Length(X[0].Value);
  SetLength(Yhat.Value, n);
  for j := 0 to n - 1 do begin
    sum := 0.0;
    for i := 0 to k - 1 do
      sum := sum + X[i].Value[j] * Betas[i];
    Yhat.Value[j] := sum;
  end;
end;

procedure ComputeGradient(const X: WMatrixType; const Y: WVectorType; const Yhat: WVectorType; const Sigma: FloatType;
  out GradBeta: CVectorType; out GradLogSigma: FloatType);
var
  i, j, n, k: Integer;
  Resid, SS: FloatType;
begin
  k := Length(X);
  n := Length(Y.Value);
  SetLength(GradBeta, k);
  for i := 0 to k - 1 do
    GradBeta[i] := 0.0;
  SS := 0.0;
  for j := 0 to n - 1 do begin
    Resid := Y.Value[j] - Yhat.Value[j];
    SS := SS + Resid * Resid;
    for i := 0 to k - 1 do
      GradBeta[i] := GradBeta[i] - X[i].Value[j] * Resid;
  end;
  for i := 0 to k - 1 do
    GradBeta[i] := GradBeta[i] / Sqr(Sigma);
  GradLogSigma := n - ss / Sqr(Sigma);
end;

procedure ConditionHessianRobust(var H: CMatrixType; p: Integer);
var
  i: Integer;
  Lambda, DiagMax: FloatType;
  Hcopy: CMatrixType;
  Okay: Boolean;
begin
  { Determine max diagonal magnitude }
  DiagMax := 0.0;
  for i := 0 to p-1 do
    if Abs(H[i, i]) > DiagMax then
      DiagMax := Abs(H[i, i]);

  if DiagMax = 0.0 then
    DiagMax := 1.0;

  Lambda := Tolerance * DiagMax;

  { Try multiple damping levels until invertible }
  Okay := False;
  while (not Okay) and (Lambda < LowTolerance) do begin
    { Copy H }
    SetLength(Hcopy, p, p);
    for i := 0 to p - 1 do
      Move(H[i, 0], Hcopy[i, 0], p * SizeOf(FloatType));

    { Add λI }
    for i := 0 to p - 1 do
      Hcopy[i, i] := Hcopy[i, i] + Lambda;

    { Test invertibility }
    Okay := TestInvertCMatrix(Hcopy);

    if not Okay then
      Lambda := Lambda * 10;  { Increase damping }
  end;

  { If ok, replace H with conditioned version }
  if Okay then
    H := Hcopy
  else begin
    { FINAL FAIL-SAFE — very strong damping }
    for i := 0 to p - 1 do
      H[i, i] := H[i, i] + Lambda;
  end;
end;

{procedure MatrixVectorMultiply(const A: CMatrixType;
                               const v: CVectorType;
                               var Result: CVectorType);
var
  i, j, n: Integer;
begin
  n := Length(v);
  SetLength(Result, n);
  for i := 0 to n - 1 do begin
    Result[i] := 0.0;
    for j := 0 to n - 1 do
      Result[i] := Result[i] + A[i,j] * v[j];
  end;
end;}

procedure UpdateNewton(const X: WMatrixType; const GradBeta: CVectorType; var Betas: CVectorType);
var
  H, InvH: CMatrixType;
  Step: CVectorType;
  i, p: Integer;
  Okay: Boolean;
begin
  p := Length(Betas);

  { Build Hessian = Xᵀ X }
  H := CreateATA(X);

  { Condition it }
  ConditionHessianRobust(H, p);

  { Invert it }
  InvH := InvertCMatrix(H, Okay);
  if not Okay then begin
    AddToErrorStack('Newton update: Hessian not invertible.');
    Exit;
  end;

  { Step = H^[-1] * grad }
  SetLength(Step, p);
  //MatrixVectorMultiply(InvH, GradBeta, Step);
  Step := MultiplyCMatrixVector(InvH, GradBeta);

  { β := β - step }
  for i := 0 to p - 1 do
    Betas[i] := Betas[i] - Step[i];
end;

procedure ComputeXtXDiag(const X: WMatrixType; out XtXDiag: CVectorType);
var
  i, j, k, n: Integer;
  d: FloatType;
begin
  k := Length(X);
  if k = 0 then Exit;

  n := Length(X[0].Value);

  // Allocate p diagonal values.
  SetLength(XtXDiag, k);

  // Initialize.
  for j := 0 to k - 1 do
    XtXDiag[j] := 0.0;

  // Accumulate X[j,i]^2 across all observations.
  for i := 0 to n - 1 do
    for j := 0 to k - 1 do begin
      d := X[j].Value[i];
      XtXDiag[j] := XtXDiag[j] + Sqr(d);
    end;
end;

procedure UpdateGD(const GradBeta: CVectorType; const XtXDiag: CVectorType;
  var Betas: CVectorType; const BaseLR: FloatType);
var
  i: Integer;
  Step: FloatType;
begin
  for i := 0 to High(Betas) do begin
    Step := BaseLR * GradBeta[i] / (XtXDiag[i] + 1e-12);

    // Clip runaway steps.
    if Step > 1.0 then Step := 1.0;
    if Step < -1.0 then Step := -1.0;

    Betas[i] := Betas[i] - Step;
  end;
end;

// Fisher Scoring is for Logistic Regression.
{procedure UpdateFisherScoring(const X: WDataMatrix;const Y: WVectorType; var Betas: CVectorType; out Step: CVectorType);
var
  n, p, i, j, k: Integer;
  eta, pi, wi: FloatType;
  W: CVectorType;         // diagonal weights p*(1-p)
  Resid: CVectorType;     // y - p
  XtW: CMatrixType;       // X'W
  XtWX: CMatrixType;      // X'WX (Hessian approx)
  InvXtWX: CMatrixType;   // inverse Hessian
  XtWResid: CVectorType;  // X'(y - p)
begin
  p := Length(IndVar);
  n := Length(IndVar[0].Value);

  // Allocate
  SetLength(W, n);
  SetLength(Resid, n);

  // Compute p_i and w_i = p(1-p)
  for i := 0 to n - 1 do begin
    // η = xᵢ·β
    eta := 0.0;
    for j := 0 to p - 1 do
      eta := eta + Betas[j] * IndVar[j].Value[i];

    // logistic p
    pi := 1.0 / (1.0 + Exp(-eta));

    // weight
    wi := pi * (1.0 - pi);
    W[i] := wi;

    // residual
    Resid[i] := DepVar.Value[i] - pi;
  end;

  // Compute XtW  (p × n matrix)
  SetLength(XtW, p);
  for j := 0 to p - 1 do begin
    SetLength(XtW[j], n);
    for i := 0 to n - 1 do
      XtW[j][i] := IndVar[j].Value[i] * W[i];
  end;

  // XtWX = XtW * X
  SetLength(XtWX, p);
  for j := 0 to p - 1 do begin
    SetLength(XtWX[j], p);
    for k := 0 to p - 1 do begin
      XtWX[j][k] := 0.0;
      for i := 0 to n - 1 do
        XtWX[j][k] := XtWX[j][k] + XtW[j][i] * IndVar[k].Value[i];
    end;
  end;

  // Invert Hessian
  CreateInverse(XtWX, InvXtWX);

  // Compute XtW * Resid
  SetLength(XtWResid, p);
  for j := 0 to p - 1 do begin
    XtWResid[j] := 0.0;
    for i := 0 to n - 1 do
      XtWResid[j] := XtWResid[j] + XtW[j][i] * Resid[i];
  end;

  // Final Fisher scoring update:
  // Step = (XtWX)^(-1) * XtWResid
  SetLength(Step, p);
  for j := 0 to p - 1 do begin
    Step[j] := 0.0;
    for i := 0 to p - 1 do
      Step[j] := Step[j] + InvXtWX[j][i] * XtWResid[i];
  end;

  // Apply update to Betas
  for j := 0 to p - 1 do
    Betas[j] := Betas[j] + Step[j];
end;}

procedure UpdateNesterov(const X: WMatrixType; const Y: WVectorType; const Betas: CVectorType;
  var NewBetas, Velocity : CVectorType; const LearningRate, MomentumCoef, Sigma: FloatType);
var
  i, j, p, n: Integer;
  TempBetas, Yhat: WVectorType;
  GradBeta: CVectorType;
  GradLogSigma: FloatType;
begin
  p := Length(Betas);
  n := Length(Y.Value);

  { --- allocate arrays --- }
  if Length(NewBetas) <> p then SetLength(NewBetas, p);
  if Length(Velocity) <> p then SetLength(Velocity, p);
  SetLength(TempBetas.Value, p);
  SetLength(GradBeta,  p);
  SetLength(Yhat.Value, n);
  Yhat.Name := 'Yhat';

  { --- 1. Look-ahead step: TempBetas = Betas + mu * Velocity --- }
  for i := 0 to p - 1 do
    TempBetas.Value[i] := Betas[i] + MomentumCoef * Velocity[i];

  { --- 2. Compute Yhat = X * TempBetas --- }
  // This is your linear predictor
 // for i := 0 to n - 1 do
   // Yhat.Value[i] := VectorDotProduct(X[i], TempBetas);
  for i := 0 to n - 1 do begin
    Yhat.Value[i] := 0.0;
    for j := 0 to p - 1 do
      Yhat.Value[i] := Yhat.Value[i] + X[j].Value[i] * TempBetas.Value[j];
end;

  { --- 3. Compute gradient at look-ahead params --- }
  ComputeGradient(X, Y, Yhat, Sigma, GradBeta, GradLogSigma);

  { --- 4. Update velocity: v = mu * v - lr * grad --- }
  for i := 0 to p - 1 do
    Velocity[i] := MomentumCoef * Velocity[i] - LearningRate * GradBeta[i];

  { --- 5. Parameter update: beta <- beta + v --- }
  for i := 0 to p - 1 do
    NewBetas[i] := Betas[i] + Velocity[i];
end;

{procedure UpdateNesterov(const X: WMatrixType; const Y: WVectorType; const Betas: CVectorType;
  var NewBetas,Velocity: CVectorType; const LearningRate, MomentumCoef: FloatType);
// Momentum Coef is μ, typically 0.8–0.95
var
  TempBetas, Grad: CVectorType;
  i, k: Integer;
begin
  k := Length(Betas);

  // Ensure correct sizes
  if Length(Velocity) <> k then begin
    SetLength(Velocity, k);
    for i := 0 to k - 1 do
      Velocity[i] := 0.0;
  end;

  SetLength(TempBetas, k);
  SetLength(Grad, k);
  SetLength(NewBetas, k);

  // Compute "look-ahead" parameters: β + μ * v
  for i := 0 to k - 1 do
    TempBetas[i] := Betas[i] + MomentumCoef * Velocity[i];

  // Evaluate the gradient at the look-ahead point
  ComputeGradient(X, Y, TempBetas, Grad);

  // Update the velocity v = μv – η∇L(β + μv)
  for i := 0 to k - 1 do
    Velocity[i] := MomentumCoef * Velocity[i] - LearningRate * Grad[i];

  // Update parameters β_new = β + v
  for i := 0 to k - 1 do
    NewBetas[i] := Betas[i] + Velocity[i];
end;}

procedure MLERegression(const X: WMatrixType; const Y: WVectorType;
  out Betas: CVectorType; out Sigma: FloatType);
var
  i, j, n, k, Iter, LocalMaxIter: Integer;
  Xmean, XStd, Velocity: CVectorType;
  Residual, YHat: WVectorType;
  GradBeta, OldBetas, XtXDiag: CVectorType;
  GradLogSigma, LL, LogSigma, OldLogSigma: FloatType;
  //Standardize: Boolean = True;
  Xs: WMatrixType;
begin
  k := Length(X);
  n := Length(Y.Value);

  // --- allocate and copy X (we will standardize optionally) ---
  SetLength(Xs, k);
  SetLength(Xmean, n);
  SetLength(XStd, n);
  BasicFuncs.DeepCopyWMatrix(X, Xs);
{  for i := 0 to k - 1 do begin        //Use DeepCopyMatrix
    Xs[i].Name := X[i].Name;
    SetLength(Xs[i].Value, n);
    for j := 0 to n - 1 do
      Xs[i].Value[j] := X[i].Value[j];
  end;}
  SetLength(Betas, k);
  SetLength(Velocity, k);
  for i := 0 to k - 1 do begin
    Betas[i] := 0.0;
    Velocity[i] := 0.0;
  end;
  if RegressionMode = MLEGD then
    for i := 0 to k - 1 do
      if Xs[i].Name <> 'Intercpt' then
        StandardizeRetainWColumn(Xs[i], Xmean[i], Xstd[i]);
  Sigma := 1.0;
  LogSigma := Ln(Sigma);
  OldLogSigma := 0.0;
  LocalMaxIter := IfThen(UserMaxIter > 0, UserMaxIter, SlowMaxIter);
  SetLength(OldBetas, k);
  SetLength(XtXDiag, k);   // Do I need?
  for  j := 0 to k - 1 do
    OldBetas[k] := 0.0;
  for Iter := 1 to LocalMaxIter do begin
    ComputePrediction(X, Betas, Yhat);
    Sigma := Exp(LogSigma);
    ComputeXtXDiag(X, XtXDiag);
    ComputeGradient(X, Y, Yhat, Sigma, GradBeta, GradLogSigma);

    case RegressionMode of
      MLENR: UpdateNewton(X, GradBeta, Betas);
      MLEGD: UpdateGD(GradBeta, XtXDiag, Betas, Sigma);
      //MLEFS: UpdateFisherScoring(X, GradBeta, Betas);
      MLENest: UpdateNesterov(X, Y, OldBetas, Betas, Velocity, 0.01, 0.9, Sigma);
    end;

    LogSigma := LogSigma - 0.001 * GradLogSigma;
    LL := LogLikelihood(Y, Yhat, Exp(LogSigma));

    if ((Iter mod 100) = 0) and (Verbose or DebugOn) then begin
      Write('Iter=', Iter : 7, '    LogL= ', LL:14:10, '   Sigma=', Sigma:13:10, '   Betas: ');
      for j := 0 to k -1 do write(Betas[j]: 8 : 10, ' ');
      writeln;
    end;
    if CheckConvergence(OldBetas, Betas, OldLogSigma, LogSigma, 1e-10) then Break;

    OldLogSigma := LogSigma;
    for i := 0 to k - 1 do
      OldBetas[i] := Betas[i];
  end;
  if RegressionMode = MLEGD then
    for i := 0 to k - 1 do
      if Xs[i].Name <> 'Intercpt' then
        UnstandardizeColumn(Xs[i], Xmean[i], Xstd[i]);

  // Compute Residuals.
  SetLength(Residual.Value, n);
  for i := 0 to n - 1 do
    Residual.Value[i] := Y.Value[i] - YHat.Value[i];

  // Iteration Halting Report.
  Sigma := Exp(LogSigma);
  if (Iter >= LocalMaxIter) then
    Writeln('Converged at maximum iteration of ', Iter)
  else
    Writeln('Converged at iteration ', Iter);
  Writeln;
  Writeln('ML Regression Results Independent Variables: ', k, ' Observations: ', n);
  if RegressionMode = MLEGD then
    writeln('Gradient Descent Method')
  else if RegressionMode = MLENR then
    writeln('Newton-Raphson Method')
  else if RegressionMode = MLENest then
    writeln('Nesterov Method');

  Writeln(StringOfChar('-', 90));
  Writeln('No. ', 'Variable' : 14, 'Beta' : 18);
  for i := 0 to k - 1 do begin
    write(i + 1 : 3, ' ');
      write(X[i].Name : 14);
    Writeln(Betas[i] : 18 : 5);
  end;
  Writeln('Sigma: ', Sigma : 18: 10);
  SaveRegressionResults(Residual, YHat);
  Pause;
end;

{ Ordinary Least Squares }
procedure MainRegression(const WVar: IVectorType);
var
  IndData: WMatrixType;
  DepVar: WVectorType;
  IndVar: CIVectorType;
  i, iDepVar, NumIndVar: Integer;
  Betas, StdErrors, PValues: CVectorType;
  Sigma, RSquared, AdjRSquared: FloatType;

procedure OLS(const X: WMatrixType; const Y: WVectorType;
  out Betas, StdErrors, PValues: CVectorType; out RSquared, AdjRSquared: Double);
var
  xp, a, b: FloatType;
  n, k, i: Integer;
  XtX: CMatrixType;
  XtXInv: CMatrixType;
  TStats, XtY: CVectorType;
  Residuals, YHat: WVectorType;
  SSE, SST, SSR, Sigma2, MeanY, SumY,
    FStat, PFStat, DFReg, DFError, MSE, MSR,
    RDW, RMean, RVariance, RStandardDeviation, RRMSE, RSkewness, RKurtosis: FloatType;
  Okay: Boolean;
begin
  // Get dimensions -- here, n rather than m is the number of rows.
  n := Length(X[0].Value);   // Number of observations.
  k := Length(X);            // Number of independent variables.

  // Initial check of dimensions.
  if (n < k) or (Length(Y.Value) <> n) then begin
    Writeln('Error on regression: invalid matrix dimensions or insufficient observations.');
    Exit;
  end;
  if Verbose or DeBugOn then
    writeln('Now starting OLS. Obs and Vars: ', n, ' ', k);

  // Initialize arrays.
  SetLength(XtX, k);
  SetLength(XtXInv, k, k);
  SetLength(Betas, k);
  SetLength(StdErrors, k);
  SetLength(PValues, k);
  SetLength(TStats, k);
  SetLength(YHat.Value, n);
  SetLength(Residuals.Value, n);
  SetLength(XtY, k);

  // Create XtX.
  XtX := CreateATA(X);
  if Verbose or DeBugOn then
    DisplayMatrixM('XtX ', XtX, T3);

  // Compute (XtX)^(-1).
  Okay := True;
  InvertCMatrixType(XtX, XtXInv, Okay);
  if Verbose or DebugOn then
    DisplayMatrixM('Inverse XtX ', XtXInv, T3);
  writeln(' xtxinv 00 ', XTXInv[0,0]);
  if not Okay then begin
    writeln('Error on Regression: matrix cannot be inverted.');
    Exit;
  end;

  // Compute XtY.
  XtY := CreateXtYC(X, Y);
  if Verbose or DeBugOn then
    DisplayVectorM('XtY ', XtY, T2);

  // Compute betas: Inv(XtX) XtY.
  Betas := MultiplyCMatrixVector(XtXInv, XtY);
  if Verbose or DeBugOn then
    DisplayVectorM('Betas ', Betas, T2);

  // Compute fitted values: YHat = Xβ.
  YHat := MultiplyWMatrixCVector(X, Betas);
  YHat.Name := 'YHat';
  FillVectorCaches(YHat);
  if Verbose or DeBugOn then
    DisplayVectorM('YHat ', YHat, T2);

  // Compute residuals: e = Y - YHat.
  for i := 0 to n - 1 do
    Residuals.Value[i] := Y.Value[i] - YHat.Value[i];
  Residuals.Name := 'Resid';
  FillVectorCaches(Residuals);
  if Verbose or DeBugOn then
    DisplayVectorM('Residuals ', Residuals, T2);

  // Compute SSE.
  SSE := 0.0;
  for i := 0 to n - 1 do
    SSE := SSE + Sqr(Residuals.Value[i]);

  // Compute mean of Y.
  SumY := 0.0;
  for i := 0 to n - 1 do
    SumY := SumY + Y.Value[i];
  MeanY := SumY / n;

  // Compute SST.
  SST := 0.0;
  for i := 0 to n - 1 do
    SST := SST + Sqr(Y.Value[i] - MeanY);

  // Compute SSR.
  SSR := SST - SSE;

  // Compute R-squared.
  if SST > 0 then
    RSquared := SSR / SST
  else
    RSquared := 0.0;

  // Compute adjusted R-squared.
  if (n > k) and (SST > 0) then
    AdjRSquared := 1 - (SSE / (n - k)) / (SST / (n - 1))
  else
    AdjRSquared := 0.0;

  // Compute Sigma2.
  Sigma2 := SSE / (n - k);

  // Compute standard errors and p-values, two-tailed.
  for i := 0 to k - 1 do begin
    if (XtXInv[i, i] * Sigma2) > 0.0 then
      StdErrors[i] := Sqrt(XtXInv[i, i] * Sigma2)
    else StdErrors[i] := 0.0;
    if StdErrors[i] > 0 then
      TStats[i] := Betas[i] / StdErrors[i]
    else
      TStats[i] := 0.0;
    PValues[i] := TDistPValue(Abs(TStats[i]), n - k);
  end;

  // Output OLS.
  Writeln;
  Writeln('OLS Regression Results Independent Variables: ', k, ' Observations: ', n);
  Writeln(StringOfChar('-', 90));
  Writeln('No. ', 'Variable' : 14, 'Beta' : 18, 'Std Error' : 18, 'T-Statistic' : 18, 'P-Value' : 18);
  for i := 0 to k - 1 do begin
    write(i + 1 : 3, ' ');
    if IndVar[i] <> -1 then
      write(WData[IndVar[i]].Name : 14)
    else
      write('Intercpt' : 14);
    Writeln(Betas[i] : 18 : Precision, StdErrors[i] : 18 : Precision, TStats[i] : 18 : Precision, PValues[i] : 18 : Precision);
  end;
  Writeln('R-squared: ', SmartFloat(RSquared));
  Writeln('Adjusted R-squared: ', SmartFloat(AdjRSquared));

  // Do ANOVA.
  DFReg := k - Ord(UseIntercept);
  DFError := n - k;
  MSR := SSR / DFReg;
  MSE := SSE / DFError;
  FStat := MSR / MSE;
    xp := (DFReg * FStat) / (DFReg * FStat + DFError);
  a := DFReg / 2;
  b := DFError / 2;
  PFStat := 1.0 - IncompleteBetaCF(xp, a, b);

  // Output ANOVA.
  Writeln;
  Writeln('Analysis of Variance');
  Writeln(StringOfChar('-', 78));
  Writeln(' ': 18, 'Sums of Squares': 24, 'Mean Squares': 24, 'DF': 12);
  Writeln('Regression          ', SmartFloat(SSR): 24, SmartFloat(MSR): 24, DFReg: 12: 0);
  Writeln('Residuals           ', SmartFloat(SSE): 24, SmartFloat(MSE): 24, DFError: 12: 0);
  Writeln('Total               ', SmartFloat(SST): 24);
  Writeln('F Statistic: ', SmartFloat(FStat));
  Writeln('Probability: ', SmartFloat(PFStat));

  // Compute residuals and associated statistics.
  RDW := DurbinWatson(Residuals);
  RMean := ArithmeticMean(Residuals);
  RRMSE := RootMeanSquare(Residuals);
  RVariance := Variance(Residuals, Sample);
  RStandardDeviation := StandardDeviation(Residuals, Sample);
  RSkewness := Skewness(Residuals, Sample);
  RKurtosis := Kurtosis(Residuals, Sample);

  // Output residuals and associated statistics.
  writeln;
  writeln('Analysis of Residuals');
  Writeln(StringOfChar('-', 66));
  Writeln('Durbin-Walson Statistic: ', SmartFloat(RDW));
  Writeln('Mean of Residuals: ', SmartFloat(RMean));
  Writeln('Root Mean Squared Error of Residuals: ', SmartFloat(RRMSE));
  Writeln('Variance of Residuals: ', SmartFloat(RVariance));
  Writeln('Standard Deviation: ', SmartFloat(RStandardDeviation));
  Writeln('Skewness of Residuals: ', SmartFloat(RSkewness));
  Writeln('Kurtosis of Residuals: ', SmartFloat(RKurtosis));

  SaveRegressionResults(Residuals, YHat);

  // Additional information.
  Writeln(TwoTailed);
end;

procedure AppendMatrixOnes(var A: WMatrixType);
var
  j, m, n: Integer;
begin
  n := Length(A);
  m := Length(A[0].Value);

  // New length of A will be m, n + 1.
  SetLength(A, n + 1);

  // Add column of 1's.
  SetLength(A[n].Value, m);
  for j := 0 to m - 1 do
    A[n].Value[j] := 1.0;
  // Name is set in calling proc.
end;

// MainRegression procedure.
begin
  // Get the dependent variable.
  iDepVar := WVar[0];

  // Get the independent variables, which start at WVar[1]. IndVar starts at zero.
  NumIndVar := Length(WVar) - 1;
  SetLength(IndVar, NumIndVar);
  for i := 0 to NumIndVar do
    IndVar[i] := WVar[i + 1];

  // Create IndData, which will house the independent variables.
  SetLength(IndData, NumIndVar);

  // Copy the vectors of WData into IndData.
  for i := 0 to NumIndVar - 1 do begin
    SetLength(IndData[i].Value, Length(WData[IndVar[i]].Value));
    DeepCopyVector(WData[IndVar[i]], IndData[i]);
  end;

  // Add the Intercept of ones.
  if UseIntercept then begin
    AppendMatrixOnes(IndData);
    IndData[NumIndVar].Name := 'Intercpt';
    IndVar[NumIndVar] := -1;  //intercept has variable number -1
    Inc(NumIndVar);
  end;

  // Check dependent and independent variables.
  if Verbose or DebugOn then begin
    writeln('All zero-based. There are ', NumIndVar, ' independent variable(s) (including any intercept).');
    write('The independent variables are: '); for i := 0 to numIndVar - 1 do write(Indvar[i], ' ');
    writeln('The dependent variable is: ', iDepVar);
    writeln('Number of columns of IndData: ', Length(IndData));
    Pause;
  end;

  // Create DepVar, which will house the dependent variable.
  SetLength(DepVar.Value, Length(IndData[0].Value));
  DeepCopyVector(WData[iDepVar], DepVar);

  // The dependent variable cannot be in the DepVar.
  for i := 0 to NumIndVar - 1 do
    if iDepVar = IndVar[i] then begin
      writeln('Error on regression: dependent variable cannot be independendent variable.');
      Pause;
      Exit;
    end;

  // Call Estimation Procedure.
  Case RegressionMode of
    // OLS Procedure.
    MLEOLS : OLS(IndData, DepVar, Betas, StdErrors, PValues, RSquared, AdjRSquared);

    //  Gradient Descent Procedure.
    MLEGD : MLERegression(IndData, DepVar, Betas, Sigma);

    // Call Newton-Raphson Procedure.
    MLENR : MLERegression(IndData, DepVar, Betas, Sigma);

    // Call Nesterov Procedure.
    MLENest : MLERegression(IndData, DepVar, Betas, Sigma);
  end;

  DisplayErrorStack;
  if FirstPass then
    SavePartialData('Save regression to a file? (y/n) ', WVar, @MainRegression);
  Pause;
end;

end.

