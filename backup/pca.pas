unit PCA;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Correlations,
  Crt,
  DataDisplay,
  DataManager,
  Globals,
  LinearAlgebra,
  SysUtils;

procedure PCAnalysis(const WVar: IVectorType);

implementation

procedure PCAnalysis(const WVar: IVectorType);
var
  i, m, n: Integer;
  CentData: WMatrixType;
  EigenVectors, CovMatrix, CorrMatrix: CMatrixType;
  EigenValues: CVectorType;
  nVar: Integer;
  MinEig, SpecRad: FloatType;
begin
  ClrScr;
  writeln('Principal Component Analysis');
  m := Length(WData[0].Value);  // This is duped in advanced section.
  n := Length(WData);
  nVar := Length(WVar);
  writeln;
  writeln('There are ', n, ' variables and ', m,' observations.');

  // Need to create new matrix, CentData.
  SetLength(CentData, nVar);
  for i := 0 to nVar - 1 do begin
    SetLength(CentData[i].Value, Length(WData[WVar[i]].Value));
    DeepCopyVector(WData[WVar[i]], CentData[i]);
  end;

  if DebugOn then
    DisplayMatrixM('CentData', CentData, T2);
  Pause;

  // Create eigenvectors and eigenvalues.
  SetLength(EigenVectors, nVar);
  SetLength(EigenValues, nVar);

  // Perform PCA with correlation or covariance matrix.
  case PCAOption of
    CorrPCA: begin
      // Center vectors of CentData.
      CenterMatrix(CentData);

      // Create covariance matrix.
      SetLength(CovMatrix, nVar);
      CreateCovarianceMatrix(CentData, CovMatrix);
      if DebugOn then
        DisplayMatrixM('Covariance', CovMatrix, T2);

      // Create eigenvectors and eigenvalues from Jacobi.
      JacobiEigenDecomposition(CovMatrix, EigenVectors, EigenValues, SpecRad);
    end;

    CovPCA: begin
      // Create correlation matrix.
      SetLength(CorrMatrix, nVar);
      CreateCorrelationMatrix(CentData, CorrMatrix);
      if DebugOn then
        DisplayMatrixM('Correlation', CorrMatrix, T2);

      // Create eigenvectors and eigenvalues from Jacobi.
      JacobiEigenDecomposition(CorrMatrix, EigenVectors, EigenValues, SpecRad);
    end;
  end;

  // Report results.
  ReportPSD(EigenValues, MinEig);

  // Save results.
  SavePartialData('Save PCA to a file? (y/n) ', WVar, @PCAnalysis);
  Pause;

end;
end.

