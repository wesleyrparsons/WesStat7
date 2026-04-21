unit Globals;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
    Crt,
    DateUtils,
    Math,
    SysUtils;

// === Precision Selection ===
// Uncomment ONE of the following or define via compiler command line:
// {$DEFINE USE_SINGLE}
  {$DEFINE USE_DOUBLE}
// {$DEFINE USE_EXTENDED}

  {$IFDEF USE_SINGLE}
  type
    FloatType = Single;
    DisplayFormatType = (T1, T2, T3);
  const
    MachineEpsilon = 1.1920928955078125e-7;
    CompilerEpsilon = 1e-6;
    CompilerHighTolerance = 1e-6;
    CompilerTolerance = 1e-6;
    CompilerLowTolerance = 1e-4;
    MinRealNumber = 1.40129846e-45;              // Not used.
    MaxRealNumber = 3.40282347e+38;              // Not used.
    VeryTinyNumber  = 1e-40;                     // Not used.
  var DefaultDisplayFormat: DisplayFormatType = T1;
    {$ENDIF}

  {$IFDEF USE_DOUBLE}
  type
    FloatType = Double;
    DisplayFormatType = (T1, T2, T3);
  const
    MachineEpsilon = 2.225e-16;
    CompilerEpsilon = 1e-12;
    CompilerHighTolerance = 1e-12;
    CompilerTolerance = 1e-9;
    CompilerLowTolerance = 1e-6;
    MinRealNumber  = 2.2250738585072020E-308;    // Not used.
    MaxRealNumber  = 1.7976931348623150E+308;    // Not used.
    VeryTinyNumber  = 1e-300;                    // Not used.
  var DefaultDisplayFormat: DisplayFormatType = T2;
  {$ENDIF}

  {$IFDEF USE_EXTENDED}
  type
    FloatType = Extended;
    DisplayFormatType = (T1, T2, T3);
    DisplayFormatType = (T1, T2, T3);
  const
    ExtendedMachineEpsilon = 1.084202172485504e-19;
    CompilerEpsilon = 1e-18;
    CompilerHighTolerance = 1e-18;
    CompilerTolerance = 1e-16;
    CompilerLowTolerance = 1e-12;
    MinRealNumber = 3.6451995318824746e-4951;    // Not used.
    MaxRealNumber = 1.1897314953572318e+4932;    // Not used.
    VeryTinyNumber  = 1e-4900;                   // Not used.
  var DefaultDisplayFormat: DisplayFormatType = T3;
  {$ENDIF}

type
  IVectorType = Array of Integer;
  CVectorType = Array of FloatType;
  CIVectorType = Array of Integer;
  SVectorType = Array of String;
  BVectorType = Array of Boolean;
  CMatrixType = Array of Array of FloatType;
  TMessageType = Array[0..99] of String;
  MissingDataType = (MDMean, MDMedian, MDZero);
  RegressionType = (MLEOLS, MLEGD, MLENR, MLEFS, MLENest);
  CorrType = (Pearson, Spearman, Kendall, Hoeffding, MI, NMI, Covar);
  PCAType = (CorrPCA, CovPCA);
  EntropyType = (SqRoot, Rice, Sturges, FD, Scott, Fixed);
  WVectorType = record
    Name: String;
    Median, Mean, StdDev, Min, Max: FloatType;
    Value: CVectorType;           // Array of FloatType.
    IsNan: BVectorType;           // Array of Boolean.
  end;
  WScalarType = record
    Name: String;
    Value: FloatType;             // Array of FloatType.
  end;
  type
  TokenArgType = record
    I: Integer;
    F: FloatType;
    S: String;
    V: Integer;
  end;
  WMatrixType = Array of WVectorType;
  ProcType = procedure;
  VProcType = procedure(const Col: WVectorType);
  VVProcType = procedure(const Col1, Col2: WVectorType);
  IProcType = procedure(const WVar: IVectorType);

const
  Title = 'Wes''s Statistics Version 6, November 2025';
  TwoTailed = 'Significance is two-tailed.';
  //MaxCol = 99;             // Maximum columns.
  //MaxRow = 999;            // Maximum rows.
  MaxTokens = 127;
  MaxObjects = 127;          // Max objects, 255 div 2, roughly.
  Pi       = 3.14159265358979323846264;      // pi.
  TwoPi    = 6.28318530717958647692529;      // 2*pi.
  InvSqrt2 = 0.70710678118654752440084;      // 1/sqrt(2).
  SqrtPi2  = 1.25331413731550025120788;      // sqrt(pi/2).
  Ln2      = 0.69314718055994530941723;      // natural log 2.
  Sqrt2Pi  = 2.50662827463100050241577;      // sqrt(2*pi).
  Sqrt2    = 1.41421356237309504880169;      // sqrt(2).
  Sqrt3    = 1.73205080756887729352745;      // sqrt(3).
  SD1      = 0.68268949213708609994572;      // SD 1 on normal distribution.
  SD2      = 0.95449973610364206493212;      // SD 2 on normal distribution.
  SD3      = 0.99730020393674846173921;      // SD 3 on normal distribution.
  E        = 2.71828182845904523536029;      // euler's number.
  Ln10     = 2.30258509299404568401799;      // natural log 10.
  Delta    = 4.66920160910299067185320;      // Feigenbaum constant.
  Phi      = 1.61803398874989484820459;      // golden ratio.

var
  { Errors & Debug }
  iES: Integer = 0;
  ErrorStack: TMessageType;
  DebugOn: Boolean = True;

  { Scalars }
  SData: Array of WScalarType;
  nScal: Integer = 0;

  { Main Data }
  nRow, nCol: Integer;
  WData: WMatrixType;
  WInput: String;
  FirstPass: Boolean;
  MissingData: MissingDataType = MDMean;

  { Logging }
  TeeFile: TextFile;
  TeeActive: Boolean = False;
  OldInOut: Pointer = nil;
  OldFlush: Pointer = nil;
  TeeFileName: String;

  { Objects }
  Command: String;
  WTokenArg: Array[0..MaxTokens] of TokenArgType;      // Contains Integer, Float, & String Tokens.
  WOpt: Array[0..MaxObjects] of String;
  WArg: Array[0..MaxObjects] of FloatType;
  WStr: Array[0..MaxObjects] of String;
  iTA, nTA: Integer;
  nArg, nOpt, nStr: Integer;
  iArg, iOpt, iStr: Integer;
  TVN: Integer = 0;

  { Display Options }
  Verbose: Boolean = True;
  ShowSig: Boolean = True;
  LargeThreshold: Float = 1e8;
  SmallThreshold: Float = 1e-6;
  Width: Integer = 14;
  Precision: Integer = 5;
  Tab: Integer = 50;
  ShowSignificance: Boolean = False;
  FilePath: String;
  Logging: Boolean = False;
  LogFile: Text;

  { Statistics Options }
  UseIntercept: Boolean = True;
  Ext: Boolean = True;
  Sample: Boolean = True;
  PCAOption: PCAType = CorrPCA;
  Corr: CorrType;
  RegressionMode: RegressionType = MLEGD;
  Bins: EntropyType = SqRoot;
  Reps: Integer = 1000;
  SampleSize: Integer = 30;
  SampleCutoff: Integer = 2000;
  Permutations: Integer = 500;
  BinCount: Integer = 10;

  { Iterations & Tolerance }
  HighTolerance, Tolerance, LowTolerance, UserTolerance: FloatType;
  UserMaxIter: Integer;
  FastMaxIter: Integer = 50;           // Newton, refinement.
  MediumMaxIter: Integer = 400;        // LM, QR, eigen.
  SlowMaxIter: Integer = 2000;         // Gradient / EM.
  VerySlowMaxIter: Integer = 5000;

procedure DisplayErrorStack;
procedure DisplayErrorStackM;
procedure AddToErrorStack(const Mess: String);
procedure SuppressErrorStack;
procedure CheckColumn(const Mess: String; const Col: WVectorType); overload;
procedure CheckColumn(const Mess: String; const Col: CVectorType); overload;
procedure ReportProgramState;

implementation

procedure Pause;
begin
  Write('Hit <CR> to continue.');
  readln;
end;

procedure DisplayErrorStack;
var
  k: Integer;
begin
  for k := 0 to iES - 1 do
    writeln(ErrorStack[k]);
  iES := 0;
end;

procedure DisplayErrorStackM;
begin
  write('There are ', iES, ' errors.');
  if iES > 0 then begin
    writeln('They are: ');
    DisplayErrorStack;
  end;
  writeln;
end;

procedure AddToErrorStack(const Mess: String);
var
  k: Integer;
begin
  if iES > 99 then begin
    WriteLn('Errors in stack exceed 99; further errors will not be added.');
    Exit;
  end
  else if iES > 0 then begin
    for k := 0 to iES - 1 do
      if Mess = ErrorStack[k] then Exit;
    ErrorStack[iES] := Mess;
    Inc(iES);
  end;
end;

procedure SuppressErrorStack;
begin
  iES := 0;
end;

{ Display Columns -- depends on no other units }
procedure CheckColumn(const Mess: String; const Col: WVectorType); overload;
var
  j, m: Integer;
begin
  if WData = nil then Exit;
  m := Length(Col.Value);
  writeln(Mess, '  Column:    ', Col.Name, '   Rows: ', m);
  for j := 0 to m - 1 do
    write('  ', Col.Value[j]: 5: 2);
  writeln;
  write('Mean = ', Col.Mean : 5 : 2, '  ');
  write('Median = ', Col.Median : 5 : 2, '  ');
  write('Std. Dev. = ', Col.StdDev : 5 : 2, '  ');
  write('Minimum = ', Col.Min : 5 : 2, '  ');
  write('Maximum = ', Col.Max : 5 : 2);
  writeln;
end;

procedure CheckColumn(const Mess: String; const Col: CVectorType); overload;
var
  j, m: Integer;
begin
  if WData = nil then Exit;
  m := Length(Col);
  writeln(Mess, '   Rows: ', m);
  for j := 0 to m - 1 do
    write('  ', Col[j]: 5: 2);
  writeln;
end;

procedure ReportProgramState;
var
  StartTime: TDateTime;
  i: Integer;
  MemUsed: Cardinal;
  EnvStr: String;
  Locale: String;
  BuildInfo: String;
begin
  ClrScr;
  StartTime := Now;

  writeln(Title);
  writeln('Source code by Wesley R. Parsons,');
  writeln('wespar@bellouth.net, www.wespar.com, Miami, Florida, 2025.');


  // Build and environment details
  BuildInfo := 'Free Pascal ' + {$I %FPCVERSION%} + ' on ' + {$I %FPCTARGETOS%} + '/' + {$I %FPCTARGETCPU%};
  writeln('Build Info: ', BuildInfo);
  writeln('Executable: ', ParamStr(0));
  writeln('Current Directory: ', GetCurrentDir);
  writeln('Command Line Args: ', ParamCount);

  // Runtime and resource utilization
  MemUsed := GetHeapStatus.TotalAllocated;
  writeln('Memory Used: ', MemUsed, ' bytes');
  writeln('Uptime: ', FormatDateTime('hh:nn:ss', Now - StartTime));

  // Configuration and localization
  Locale := GetEnvironmentVariable('LANG');
  if Locale = '' then Locale := 'Unknown';
  writeln('Locale: ', Locale);
  EnvStr := GetEnvironmentVariable('PATH');
  writeln('PATH Length: ', Length(EnvStr));

  // Row Dump
  writeln('There are ', nRow, ' rows/observations.');

  // Variable dump
  if nCol > 0 then begin
    writeln('There are  ', nCol, ' variables/columns:');
    writeln('The variables are:');
    for i := 0 to High(WData) do begin
      write(WData[i].Name: 15);
      if (i + 1) mod 5 = 0 then writeln;
    end;
  end;
  writeln;

  // Scalar Dump
  writeln('There are ', nScal, ' scalars: ');
  for i := 0 to nScal - 1 do
    write(SData[i].Name: 15);
  writeln;

  //Errors
  writeln('There are curently ', iES, ' errors.');
  DisplayErrorStack;
end;

end.

