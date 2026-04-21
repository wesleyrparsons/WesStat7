
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

  // Run any test routines.
  TestRoutines;

  // Start in MDI or CLI.
  Facade := CLI;

  // Initialize.
  nRow := 0;
  nCol := 0;
  nScal := 0;
  WData := nil;
  SData := nil;

  // Write the introduction
  ClrScr;
  writeln;
  writeln(Title);
  write('This program requires Data. Data can come from a CSV file, ');
  writeln('as random data, as test data, or be input manually.');
  writeln;

  // Run the CLI;
  RunCLI;

{  // Start the main loop.
  while True do begin

  // Run eith CLI or MDI.
    if Facade = CLI then begin             // Run CLI.
      writeln('Command Line Interface');
      RunCLI;
    end
    else begin                             // Run MDI.
      if WData = nil then                  // First time go to ObtainData.
        repeat
          ObtainData;  // Lower limits are 1 observation and 1 variable.
          if Facade = CLI then Break;
          if (nRow < 1) or (nCol < 1) then
            writeln('Error on Data: rows or columns are < 1.');
        until (nRow > 0) and (nCol > 0);
      ConvertNaNinData;

      // Main menu.
      ClrScr;
      writeln(Title);
      writeln('Data Analysis Page');
      writeln('Menu Interface');
      ReportDataState;
      writeln;
      writeln('Data Analysis Options: ');
      writeln('  U Univariate Statistics');
      writeln('  B Bivariate Statistics and Correlations');
      writeln('  M Multivariate Statistics');
      writeln('  R Regression');
      writeln('  C Multivariate Correlations');
      writeln('  I Nominal Analaysis (Chi-Square)');
      writeln('  A Principal Component Analysis');
      writeln('  L Linear Algebra');
      writeln('  J Command Line Driven Interface');
      GeneralOptions1;
      writeln('  X Exit');
      write('>');
      readln(WInput);
      WInput := UpCase(WInput);
      case WInput of
        'U': DisplayUnivariateStatistics(1);
        'M': MultivariateStatistics;
//        'C': MultivariateCorrelations(Corr);
        'R': MainRegression;
        'I': NominalStatistics;
        'A': PCAnalysis;
        'B': BivariateCorrelations(1, 2);
        'L': MatrixOptions(WData);
        'J': Facade := CLI;
        'X': Halt;
      end;
      GeneralOptions2(WInput);
    end;}
end.

