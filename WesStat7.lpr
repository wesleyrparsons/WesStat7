
program WesStat7;

{$mode objfpc}{$H+}
{$APPTYPE CONSOLE}
{ WesStat7 September 24, 2025
  Source code by Wesley R. Parsons,
  wespar@bellouth.net, www.wespar.com, Miami, Florida, 2025.}

// Data.Value is contained in WData.Value[0..nRow-1, 0..nCol-1].
// In routines, m can be nRow-1 or High(Col) or High(A),
//   and n can be nCol-1, or High(A[0]).
// j indexes rows, and i indexes columns.
// mu is mean,  sigma is standard devision amd sigma2 is variance.
// iRow is a specific row, and iCol is a specific column.

uses
  Bivariate,
  Check,
  BasicFuncs,
  Correlations,
  Crt,
  DataDisplay,
  DataManager,
  FloatTypes,
  Globals,
  CLI,
  LinearAlgebra,
  Math,
  Multivariate,
  Nominal,
  Parser,
  PCA,
  Regression,
  SysUtils,
  Testing,
  Univariate,
  WesUnit,
  WinDisplay;

begin
  // Resize the conmsole if possible.
  ResizeConsoleBuffer;

  // Initialize.
  nRow := 0;
  nCol := 0;
  nScal := 0;
  WData := nil;
  SData := nil;

  // Write the introduction.
  ClrScr;
  writeln;
  writeln(Title);
  write('This program requires Data. Data can come from a CSV file, ');
  writeln('as random data, as test data, or be input manually.');
  writeln;

  // Run the CLI.
  RunCLI;

end.

