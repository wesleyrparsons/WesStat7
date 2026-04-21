
program WesStat7;

{$mode objfpc}{$H+}

{ WesStat7 September 24, 2025
  Source code by Wesley R. Parsons,
  wespar@bellouth.net, www.wespar.com, Miami, Florida, 2025.}

//Data.Value is contained in WData.Value[0..nRow-1, 0..nCol-1]
//In routines, m can be nRow-1 or High(Col) or High(A),
//  and n can be nCol-1, or High(A[0])
//j indexes rows, and i indexes columns
//mu is mean,  sigma is standard devision amd sigma2 is variance
//iRow is a specific row, and iCol is a specific column
//Col is the TVector being used
//WData.Value is the TMatrix being used

uses
  Bivariate,
  Check,
  BasicFuncs,
  Correlations,
  Crt,
  Data.ValueInAndOut,
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
  Univariate,
  WesUnit,
  WinDisplay;

var
  WInput: Char;
  i: Integer;
  Line: String;
procedure ObtainData.Value;
var
  WData.ValueInput: String;
  FileName: String;
begin
  ClrScr;
  writeln(Title);
  writeln('Data.Value Entry Page');             // move obtainData.Value to Data.ValueInAndOut??
  writeln;
  write('This program requires Data.Value. Data.Value can come from a CSV file, ');
  writeln('as random Data.Value, as test Data.Value, or be input manually.');
  writeln;
  writeln('Options: ');
  writeln('  C Create Random Data.Value');
  writeln('  T1 Test Data.Value for Statistics, 12 rows x 5 columns');
  writeln('  T2 Test Data.Value for Linear Algebra, 2 rows x 2 columns');
  writeln('  T3 Test Data.Value for Linear Algebra, 3 rows x 3 columns');
  writeln('  T4 Test Data.Value for Linear Algebra, 5 rows x 3 columns');
  writeln('  R Read Data.Value from CSV file');
  writeln('  M Manually Enter Data.Value');
  writeln('  O Modify Data.Value');
  writeln('  Z Check Routines');
  writeln('  X Exit');
  write('>');
  readln(WData.ValueInput);
  WData.ValueInput := UpCase(WData.ValueInput);
  case WData.ValueInput of
    'C': begin
        repeat    // max rows is 999 and max cols is 99
         write('Enter number of columns (variables): ');
         readln(nCol);
         if (nCol <= 0) or (nCol > 100) then writeln('Number is too big or too small.');
        until (nCol > 0) and (nCol < 100);
        repeat
          write('Enter number of rows (observations): ');
          readln(nRow);
          if (nCol <= 0) or (nCol > 1000) then writeln('Number is too big or too small.');
        until (nRow > 0) and (nRow < 1000);
        writeln;
        SetLength(WData.Value, nRow, nCol);
        SetLength(WData.Name, nCol);
        for i := 0 to nCol - 1 do
          WData.Name[i] := 'Var' + IntToStr(i);
        CreateRandomData.Value(WData.Value);
        Data.ValueCorrelateAndZero(WData.Value);
        writeln('Data.Value consists of ', nCol, ' Columns and ', nRow, ' Rows.');
      writeln('Data.Value created successfully.');
      Pause;
    end;
    'T1': begin
      nRow := 12;
      nCol := 5;
      SetLength(WData.Value, nRow, nCol);
      SetLength(WData.Name, nCol);
        WData.Value[0,0] := 2;   WData.Value[0,1] := 3;   WData.Value[0,2] := 4;   WData.Value[0,3] := 5;   WData.Value[0,4] := 6;
        WData.Value[1,0] := 5;   WData.Value[1,1] := 6;   WData.Value[1,2] := 9;   WData.Value[1,3] := 1;   WData.Value[1,4] := 4;
        WData.Value[2,0] := 6;   WData.Value[2,1] := 4;   WData.Value[2,2] := 3;   WData.Value[2,3] := 11;  WData.Value[2,4] := 14;
        WData.Value[3,0] := 1;   WData.Value[3,1] := 6;   WData.Value[3,2] := 5;   WData.Value[3,3] := 0;   WData.Value[3,4] := 3;
        WData.Value[4,0] := 1;   WData.Value[4,1] := 5;   WData.Value[4,2] := 3;   WData.Value[4,3] := 14;  WData.Value[4,4] := 17;
        WData.Value[5,0] := 9;   WData.Value[5,1] := 3;   WData.Value[5,2] := 8;   WData.Value[5,3] := 9;   WData.Value[5,4] := 12;
        WData.Value[6,0] := 9;   WData.Value[6,1] := 3;   WData.Value[6,2] := 16;  WData.Value[6,3] := 0;   WData.Value[6,4] := 3;
        WData.Value[7,0] := 1;   WData.Value[7,1] := 4;   WData.Value[7,2] := 9;   WData.Value[7,3] := 1;   WData.Value[7,4] := 5;
        WData.Value[8,0] := 10;  WData.Value[8,1] := 5;   WData.Value[8,2] := -3;  WData.Value[8,3] := 13;  WData.Value[8,4] := 16;
        WData.Value[9,0] := 9;   WData.Value[9,1] := 13;  WData.Value[9,2] := 18;  WData.Value[9,3] := 9;   WData.Value[9,4] := 11;
        WData.Value[10,0] := 9;  WData.Value[10,1] := 3;  WData.Value[10,2] := 16; WData.Value[10,3] := 0;  WData.Value[10,4] := 4;
        WData.Value[11,0] := 11; WData.Value[11,1] := 4;  WData.Value[11,2] := -9; WData.Value[11,3] := 1;  WData.Value[11,4] := 5;
        WData.Name[0] := 'Ht'; WData.Name[1] := 'Wt'; WData.Name[2] := 'Size'; WData.Name[3] := 'A'; WData.Name[4] := 'b';
      writeln('Data.Value created successfully.');
      Pause;
    end;
    'T2': begin
      nRow := 2;
      nCol := 2;
      SetLength(WData.Value, nRow, nCol);   //invertible
      SetLength(WData.Name, nCol);
      WData.Name[0] := 'First'; WData.Name[1] := 'Second';
      WData.Value[0,0] := 1; WData.Value[0,1] := 0;
      WData.Value[1,0] := 0; WData.Value[1,1] := 2;
      writeln('Data.Value created successfully.');
    end;
    'T3': begin
      nRow := 3;
      nCol := 3;
      SetLength(WData.Value, nRow, nCol);             //invertible det =-9
      SetLength(WData.Name, nCol);
      WData.Name[0] := 'Way'; WData.Name[1] := 'Here'; WData.Name[2] := 'Avenue';
      WData.Value[0,0] := 2; WData.Value[0,1] := 1; WData.Value[0,2] := 3;
      WData.Value[1,0] := 0; WData.Value[1,1] := 1; WData.Value[1,2] := 4;
      WData.Value[2,0] := 5; WData.Value[2,1] := 2; WData.Value[2,2] := 1;
      writeln('Data.Value created successfully.');
    end;
    'T4': begin
      nRow := 5;
      nCol := 3;
      SetLength(WData.Value, nRow, nCol);
      SetLength(WData.Name, nCol);
      WData.Name[0] := 'St'; WData.Name[1] := 'Yeah'; WData.Name[2] := 'Avenue';
      WData.Value[0,0] := 2; WData.Value[0,1] := 1; WData.Value[0,2] := 3;
      WData.Value[1,0] := 0; WData.Value[1,1] := 1; WData.Value[1,2] := 4;
      WData.Value[2,0] := 9; WData.Value[2,1] := 2; WData.Value[2,2] := 1;
      WData.Value[3,0] := 2; WData.Value[3,1] := 7; WData.Value[3,2] := 2;
      WData.Value[4,0] := 1; WData.Value[4,1] := 6; WData.Value[4,2] := 1;
      writeln('Data.Value created successfully.');
    end;
    'R': begin
      write('Enter CSV filename (for example, matrix.csv): ');
      readln(FileName);
      ReadCSVData.Value(FileName, WData.Value, nRow, nCol);
      writeln('Read ', FileName, ' with ', nRow, ' rows and ', nCol, ' columns.');
    end;
    'M': begin
      ManuallyInputData.Value(WData.Value, nRow, nCol);
    end;
    'O': begin
      ModifyData.Value(WData.Value);
      nRow := Length(WData.Value);
      nCol := Length(WData.Value[0]);
    end;
    'Z': //Check routines go here
      CheckDeleteRow(WData.Value);
    'X': Exit;
  end;
  //Put variable numbers in the Labels for the working variables
  SetLength(WVar, Length(WData.Value));
  for i := 0 to High(WData.Value) do WVar[i] := i + 1;
  nVar := nCol;
end;

begin
  ResizeConsoleBuffer;
  while True do begin
    ClrScr;
    repeat
      ObtainData.Value;
      if (nRow < 2) or (nCol < 1) then
        writeln('Too little Data.Value to proceed.');
    until (nRow > 1) and (nCol > 0);
    repeat
{      Testexamples;
      Pause;
      repeat
        write('>');
        readln(Line);
        EvalLine(Line);
      until Line ='x';}
      ClrScr;
      writeln(Title);
      writeln('Data.Value Analysis Page');
      writeln;
      writeln('Data.Value Analysis Options: ');
      writeln('  U Univariate Statistics');
      writeln('  B Bivariate Statistics and Correlations');
      writeln('  M Multivariate Statistics');
      writeln('  R Regression');
      writeln('  C Multivariate Correlations');
      writeln('  I Chi-Square Analysis');
      writeln('  A Principal Component Analysis');
      writeln('  L Linear Algebra');
      writeln('  P Display Data.Value');
      writeln('  E Display Data.Value Condensed');
      writeln('  O Obtain or Modify Data.Value');
      writeln('  S Save Data.Value');
      writeln('  H Help');
      writeln('  X Exit');
      write('>');
      readln(WInput);
      WInput := UpCase(WInput);
      case WInput of
        'U': UnivariateStatistics(WVar);
        'M': MultivariateStatistics(WData.Value);
        'C': MultivariateCorrelations(WData.Value);
        'I': NominalStatistics(WData.Value);
        'A': PCAnalysis(WData.Value);
        'B': BivariateCorrelations(WVar);
        'L': MatrixOptions(WData.Value);
        'P': DisplayData.ValueA(WData.Value, WData.Name);//DisplaySelectedData.Value(WData.Value, False);
        'E': DisplaySelectedData.Value(WData.Value, True);
        'R': OLSRegression(WData.Value);
        'O': ObtainData.Value;
        'S': SaveData.Value('main', WData.Value);
        'H': writeln('No help available.');
        'X': halt;
      end;
    until WInput = 'X';
  end;
end.

