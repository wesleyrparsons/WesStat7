unit DataManager;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Crt,
  Globals,
  Math,
  SysUtils,
  WesUnit;

{ Input Data }
function IsProperVariableName(const TestName: String): Boolean;
function IsProperScalarName(const TestName: String): Boolean;
function IsProperName(const TestName: String): Boolean;
procedure ManuallyInputData;
{procedure GeneralOptions1;
procedure GeneralOptions2(const WInput: String);
function ReadPositiveInteger(const Pr: String; MaxInt: Integer): Integer;
function CheckRangedNumber(const S: String; const x: Integer; const R1, R2: Integer): Boolean; overload;
function CheckRangedNumber(const S: String; const x: FloatType; const R1, R2: FloatType): Boolean; overload;
procedure ReadVariableNumber(var iCol: Integer);
procedure ReadRowNumber(var iRow: Integer);
procedure ReadVariableNameOrNumber(var iCol: Integer);
procedure SelectColumns(const Data: WMatrixType; var WVar: IVectorType; var nVar: Integer);}

{ Create Data }
procedure CreateRandomData;
procedure DataCorrelateAndZero;
procedure GetTestData(const Str: String);
procedure FillVectorCaches(var Col: WVectorType);
procedure FillWDataCaches;

{ Access Data }
procedure PutWVector(const Col: WVectorType; const iCol: Integer);
procedure GetWVector(var Col: WVectorType; const iCol: Integer);
procedure AddWVector(const NewName: string; const NewVar: WVectorType; var AddError: Boolean);
procedure AddWConstant(const NewName: string; const NewVar: WVectorType; var AddError: Boolean);
procedure FillScalars;

{ Save & Read Data }
procedure SavePartialData(const Mess: String; PProc: ProcType); overload;
procedure SavePartialData(const Mess: String; const Col: WVectorType; VProc: VProcType); overload;
procedure SavePartialData(const Mess: String; const Col1, Col2: WVectorType; VVProc: VVProcType); overload;
procedure SavePartialData(const Mess: String; const WVar: IVectorType; IProc: IProcType);  overload;
procedure SaveData(const FileName: String);
procedure SaveFile(const X: CMatrixType; const FileName: String); overload;
procedure SaveFile(const X: WMatrixType; const FileName: String); overload;
procedure ReadCSVData(const FileName: String; var RowCount, ColCount: Integer);
procedure StartTeeOutput(const FileName: String);
procedure StopTeeOutput;
//procedure LogToFile(StartLogging: Boolean; const FileName: string = '');
//procedure LogData(Enable: Boolean; const FileName: string = '');

{ Modify Data }
procedure DeleteVariable(const iCol: Integer);
procedure DeleteScalar(const iScal: Integer);  // Zero based
procedure DeleteRow(const iRow: Integer);
procedure DataTrim(const iCol: Integer; const L, R: FloatType);
procedure WinsorizeWColumn(var Col: WVectorType; LowerPct, UpperPct: FloatType);
procedure YeoJohnsonTransform(var Data: WVectorType; Lambda: FloatType);

{Centering and Standardizing Routines}
procedure CenterColumn(var Col: WVectorType); overload;
procedure CenterColumn(var Col: CVectorType); overload;
procedure StandardizeColumn(var Col: WVectorType); overload;
procedure CenterMatrix(const Data: WMatrixType); overload;
procedure StandardizeMatrix(const Data: WMatrixType); overload;
procedure StandardizeMatrix(const Data: CMatrixType); overload;
procedure StandardizeRetainWColumn(var Col: WVectorType; var mu, Sigma: FloatType);
procedure UnstandardizeColumn(var Col: WVectorType; const mu, Sigma: FloatType);

{ Main Data Procedure }
//procedure ObtainData;

implementation

{function ReadPositiveInteger(const Pr: String; MaxInt: Integer): Integer;
var
  InputStr: String;
  Value, Code: Integer;
begin
  repeat
    Write(Pr);
    ReadLn(InputStr);

    // try to convert string to integer
    Val(InputStr, Value, Code);
    if (Code <> 0) or (Value <= 0) then
      Writeln('Please enter a valid whole number greater than 0.');
    if Value > MaxInt then begin
      Writeln('Please enter a number equal to or less than ', MaxInt, '.');
      Code := -99;
    end;
  until (Code = 0) and (Value > 0);

  Result := Value;
end;

function CheckRangedNumber(const S: String; const x: Integer; const R1, R2: Integer): Boolean; overload;
begin
  Result := False;
  if (x >= R1) and (x <= R2) then
    Result := True
  else
    writeln(S, ' is not in the required range ', R1 : 3, ' to ', R2 : 3, '.');
end;

function CheckRangedNumber(const S: String; const x: FloatType; const R1, R2: FloatType): Boolean; overload;
begin
  Result := False;
  if (x >= R1) and (x <= R2) then
    Result := True
else
  writeln(S, ' is not in the required range ', R1 : 3, ' to ', R2 : 3, '.');
end;}

// Check whether potential name is not is not in use as scalar.
function IsProperVariableName(const TestName: String): Boolean;
var
  i: Integer;
begin
  // Check if it is well-formed.
  Result := IsProperName(TestName);

  // Check not scalar name (for variables).
  for i := 0 to nScal - 1 do
    if UpCase(TestName) = UpCase(SData[i].Name) then begin
      write('Name already in use.');
      Exit(False);
    end;
end;

// Check whether potential name is not is not in use as variable.
function IsProperScalarName(const TestName: String): Boolean;
var
  i: Integer;
begin
  // Check if it is well-formed.
  Result := IsProperName(TestName);

  // Check not variable name (for scalara).
  for i := 0 to nCol - 1 do
    if UpCase(TestName) = UpCase(WData[i].Name) then begin
      write('Name already in use.');
      Exit(False);
    end;
end;

// Check whether potential name is well-formed.
function IsProperName(const TestName: String): Boolean;
  var
    i: Integer;
  begin
    // Check if first letter is alpha.
    if not (TestName[1] in ['A'..'Z', 'a'..'z']) then begin
      write('Name must begin with letter. ');
      Exit(False);
    end;

  // Check if all lettets are alphanumeric plus _.
  for i := 1 to Length(TestName) do
    if not (TestName[i] in ['A'..'Z', 'a'..'z', '_']) then begin
      write('Character ', TestName[i], ' not allowed in name. ');
      Exit(False);
    end;

  // Check not "VAR", which is reserved.
  if UpCase(Copy(TestName, 1, 3)) = 'VAR' then begin
    write('Var not allowed to begin name. ');
    Exit(False);
  end;

  // Check it is not reserved from Parserr, such as abs, or ln.

  Result := True;
end;

{procedure ReadVariableNameOrNumber(var iCol: Integer);
var                               // iCol is zero based
  InputStr: String;               // But user var numbers are one based
  i, Value, Code: Integer;        // This converts to zero based
  Found: Boolean;
begin
  Found := False;
  repeat
    write('Input variable name or number: ');
    ReadLn(InputStr);
    Val(UpCase(InputStr), Value, Code);  // try to convert string to integer
    if Code = 0 then begin       // if it is a number
//      writeln('input ', Value, Code);
      if (Value > 0) and (Value <= nCol) then begin  // and it's a variable
        iCol := Value - 1;       // to correspond with Pascal 0 array start
        Found := True;
      end
      else writeln('Variable number not found.');
    end
    else begin                   // check if it's a variable name
      for i := 0 to nCol - 1 do
        if UpCase(WData[i].Name) = UpCase(InputStr) then begin
          iCol := i;         // to correspond with Pascal 0 array start
          Found := True;
          Break;
        end;
    end;
    if not Found then writeln('Variable name not found.');
  until Found;
end;

procedure ReadScalarNameOrNumber(var iScal: Integer);
var                               // iCol is zero based
  InputStr: String;               // But user var numbers are one based
  i, Value, Code: Integer;
  Found: Boolean;
begin
  Found := False;
  repeat
    write('Input scalar name or number: ');
    ReadLn(InputStr);
    Val(UpCase(InputStr), Value, Code);  // try to convert string to integer
    if Code = 0 then begin       // if it is a number
//      writeln('input ', Value, Code);
      if (Value > 0) and (Value <= nScal) then begin  // and it's a variable
        iScal := Value - 1;       // to correspond with Pascal 0 array start
        Found := True;
      end
      else writeln('Constant number not found.');
    end
    else begin                   // check if it's a scalar name
      for i := 0 to nScal - 1 do
        if UpCase(SData[i].Name) = UpCase(InputStr) then begin
          iScal := i;         // to correspond with Pascal 0 array start
          Found := True;
          Break;
        end;
    end;
    if not Found then writeln('Scalar name not found.');
  until Found;
end;

procedure ReadVariableNumber(var iCol: Integer);
begin
  repeat
    write('Input number of variable: ');
    readln(iCol);      // iCol is the variable of interest
  until not((iCol < 1) or (iCol > nCol));
  Dec(iCol);         // to correspond with Pascal 0 array start
end;

procedure ReadRowNumber(var iRow: Integer);
begin
  repeat
    write('Input number of row/observation: ');
    readln(iRow);      // iRow is the variable of interest
  until not((iRow < 1) or (iRow > nRow));
  Dec(iRow);         // to correspond with Pascal 0 array start
end;

procedure SelectColumns(const Data: WMatrixType; var WVar: IVectorType; var nVar: Integer);
var
  i, k, n: Integer;
  VInput: Integer;
  Str: String;
begin
  n := Length(Data);
  nVar := 0;
  writeln('Enter variable numbers. Hit Enter to end. Enter ''a'' for all.');
  repeat

    //Get input
    Write('>');
    Readln(Str); //I need the ln on read
    Val(Str, VInput, k);

    //Exit on <CR>; all vars have been entered
    if Str = '' then
      if nVar = 0 then
        writeln('Input required.')
      else Break;

    //If 'A' then load all vars and exit
    if UpCase(Str) = 'A' then begin
      nVar := n;  //nVar is number of variables, up to n - 1
      SetLength(WVar, nVar);
      for k := 0 to nVar - 1 do //the Vars are 1..nVar
        WVar[k] := k + 1;
      Break;
    end;

    //if Letter other than A then continue
    if (k <> 0) and (Str <> '') then begin  //merge wioth str =''?
      WriteLn('Error: Invalid input.');
      Continue; //start loop again;
    end;

    //Check for range 1..n, if not continue
    if (VInput < 1) or (VInput > n) then begin
      WriteLn(Vinput, 'Error: Enter a variable number between 1 and ', n, '.');
      Continue; //start loop again;
    end;

    // Check if already already selected
    if nVar >= 1 then
      for k := 1 to nVar do
        if WVar[k - 1] = VInput then begin
          WriteLn('Error: Variable ', k, ' has already been selected.');
          Continue; //start loop again
        end;

    // All Ok
    SetLength(WVar, nVar + 1);
    WVar[nVar] := VInput;
    Inc(nVar);
  until (nVar = n);  //loop until load all vars
  //writeln;
  for i := 0 to nVar - 1 do Dec(WVar[i]);  //WVar should start at 0

  // Check variables selected
  if DebugOn then begin
    write('Selected ', nVar, ' variables: ');
    for i := 0 to nVar - 1 do
      write(' ', WVar[i]);
    writeln;
  end;
end;}

procedure ManuallyInputData;
var
  i, j: Integer;
begin
  write('Number of columns (variables): ');
  readln(nCol);
  write('Number of rows (observations): ');
  readln(nRow);
  SetLength(WData, nCol);
  writeln('Enter each piece of Data, follow by <CR>.');
  for i := 0 to nCol do
    SetLength(WData[i].Value, nRow);
  for j := 0 to nRow - 1 do begin
    write('Row ', j + 1);
    for i := 0 to nCol - 1 do begin
      write(' Col ', i + 1, ': ');
      read(WData[i].Value[j]);
//      write(' ');
    end;
//    writeln;
  end;
  for i := 0 to nCol - 1 do begin
    write('Enter name of variable ', i, ': ');
    ReadLn(WData[i].Name);
  end;
  // Deal with NaN.
  ConvertNaNinData;
  FillWDataCaches;
  writeln(' There are ', nRow, ' rows and ', nCol , nCol, ' columns.');
  writeln('Data were successfully entered.');
  Pause;
end;

{ Create Data }
// Use Windows function to seed random number generator.
function RtlGenRandom(pBuffer: Pointer; dwLength: Cardinal): Boolean;
  stdcall; external 'advapi32.dll' name 'SystemFunction036';

// Random number generator.
function SecureRandomFloat: FloatType;
var
  raw: Cardinal;
begin
  if not RtlGenRandom(@raw, SizeOf(raw)) then
    AddToErrorStack('Error or secure random generation.');
  Result := raw / High(Cardinal);  // Normalize to [0,1)
end;


// Create one random normal float with mean and SD.
function SecureRandomNormal(const Mean, SD: FloatType): FloatType;
var
  u1, u2, z: FloatType;
begin
  repeat
    u1 := SecureRandomFloat;
  until u1 > 0.0;
  u2 := SecureRandomFloat;
  z := Sqrt(-2 * Ln(u1)) * Cos(2 * Pi * u2);
  Result := Mean + z * SD;
end;

// In test data, correlate columns 4 and 5, and make columns 2 and 3 > 0.
procedure DataCorrelateAndZero;
var
  j: Integer;
begin
  if nCol > 4 then
    for j := 0 to High(WData[1].Value) do begin
      WData[1].Value[j] := abs(WData[1].Value[j]) + 8.0;  //  col 2 > 0
      WData[2].Value[j] := abs(WData[1].Value[j]) + 3.3;  //  col 3 > 0
      WData[4].Value[j] := WData[3].Value[j] + Random;    //cols 4 & 5 correlated
    end;
end;

procedure CreateRandomData;
var
  i, j: Integer;
begin
  try
    SetLength(WData, nCol);
    for i := 0 to nCol - 1 do begin
      WData[i].Name := 'Var' + IntToStr(i + 1);
      SetLength(WData[i].Value, nRow);
    end;
    for i := 0 to High(WData) do
      for j := 0 to High(WData[i].Value) do
        WData[i].Value[j] := SecureRandomNormal(10, 5);
    DataCorrelateAndZero;
    FillWDataCaches;
    writeln('Data consists of ', nCol, ' Columns and ', nRow, ' Rows.');
    writeln('Data created successfully, any variables Var4 & Var5 correlated, Var2 & Var3 > 0.');
  except
    writeln('Data creation unsuccessful. Reduce number of variables or observations.');
  end;
end;

// Add test data of various kinds.
procedure GetTestData(const Str: String);
var
  i: Integer;
begin
  Case UpCase(Str) of
    'T1': begin
      nRow := 12;
      nCol := 5;
      SetLength(WData, nCol);
      for i := 0 to nCol - 1 do
        SetLength(WData[i].Value, nRow);
      WData[0].Value[0] := 2;   WData[0].Value[1] := 5;   WData[0].Value[2] := 6;   WData[0].Value[3] := 1;   WData[0].Value[4] := 1;
      WData[0].Value[5] := 9;   WData[0].Value[6] := 9;   WData[0].Value[7] := 1;   WData[0].Value[8] := 10;  WData[0].Value[9] := 9;
      WData[0].Value[10] := 9;  WData[0].Value[11] := 11;

      WData[1].Value[0] := NaN;   WData[1].Value[1] := 6.6;   WData[1].Value[2] := 4;   WData[1].Value[3] := 6;   WData[1].Value[4] := 5;
      WData[1].Value[5] := 3;   WData[1].Value[6] := 3;   WData[1].Value[7] := 4;   WData[1].Value[8] := 5;   WData[1].Value[9] := 13;
      WData[1].Value[10] := 3;  WData[1].Value[11] := 4;

      WData[2].Value[0] := 4;   WData[2].Value[1] := NaN;   WData[2].Value[2] := 3;   WData[2].Value[3] := 5;   WData[2].Value[4] := 3;
      WData[2].Value[5] := 8;   WData[2].Value[6] := 16;  WData[2].Value[7] := 9;   WData[2].Value[8] := -3;  WData[2].Value[9] := 18;
      WData[2].Value[10] := 16; WData[2].Value[11] := -9;

      WData[3].Value[0] := 5;   WData[3].Value[1] := 11.111111111;   WData[3].Value[2] := 11;  WData[3].Value[3] := 0;   WData[3].Value[4] := 14;
      WData[3].Value[5] := 9;   WData[3].Value[6] := 0;   WData[3].Value[7] := 1;   WData[3].Value[8] := 13;  WData[3].Value[9] := 9;
      WData[3].Value[10] := 0;  WData[3].Value[11] := 1;

      WData[4].Value[0] := 6.1;   WData[4].Value[1] := 4.44;   WData[4].Value[2] := 14;  WData[4].Value[3] := 3;   WData[4].Value[4] := 17;
      WData[4].Value[5] := 12;  WData[4].Value[6] := 3;   WData[4].Value[7] := 5;   WData[4].Value[8] := 16;  WData[4].Value[9] := 11;
      WData[4].Value[10] := 4;  WData[4].Value[11] := 5;
      WData[0].Name := 'Ht'; WData[1].Name := 'Wt'; WData[2].Name := 'Size'; WData[3].Name := 'A'; WData[4].Name := 'baxter56';
      writeln('Data created successfully. 12 x 5. Vars 4 and 5 correlated. Two data are NaN.');
    end;
    'T9': begin
      nRow := 12;
      nCol := 5;
      SetLength(WData, nCol);
      for i := 0 to nCol - 1 do
        SetLength(WData[i].Value, nRow);
      WData[0].Value[0] := 557682;   WData[0].Value[1] := 5;   WData[0].Value[2] := 6;   WData[0].Value[3] := 1;   WData[0].Value[4] := 1;
      WData[0].Value[5] := 444449;   WData[0].Value[6] := 9;   WData[0].Value[7] := 1;   WData[0].Value[8] := 10;  WData[0].Value[9] := 9;
      WData[0].Value[10] := 9;  WData[0].Value[11] := 14441;

      WData[1].Value[0] := NaN;   WData[1].Value[1] := 6.6;   WData[1].Value[2] := 4;   WData[1].Value[3] := 6;   WData[1].Value[4] := 5;
      WData[1].Value[5] := 3456;   WData[1].Value[6] := 3;   WData[1].Value[7] := -48765487658764;   WData[1].Value[8] := 5;   WData[1].Value[9] := 13;
      WData[1].Value[10] := 3;  WData[1].Value[11] := 4444;

      WData[2].Value[0] := 4;   WData[2].Value[1] := NaN;   WData[2].Value[2] := 3.1415926;   WData[2].Value[3] := 2.7182818281115;   WData[2].Value[4] := 3;
      WData[2].Value[5] := -3763876;   WData[2].Value[6] := 16;  WData[2].Value[7] := 9.58759875;   WData[2].Value[8] := -9993;  WData[2].Value[9] := 18;
      WData[2].Value[10] := 16; WData[2].Value[11] := -9;

      WData[3].Value[0] := 2225;   WData[3].Value[1] := 11;   WData[3].Value[2] := 11;  WData[3].Value[3] := 0;   WData[3].Value[4] := 14;
      WData[3].Value[5] := 9;   WData[3].Value[6] := 0;   WData[3].Value[7] := 1.57965;   WData[3].Value[8] := 13000;  WData[3].Value[9] := 9;
      WData[3].Value[10] := 0;  WData[3].Value[11] := 1;

      WData[4].Value[0] := 644.109498074987;   WData[4].Value[1] := 8734687487648764.44;   WData[4].Value[2] := 14;  WData[4].Value[3] := 978463;   WData[4].Value[4] := 17;
      WData[4].Value[5] := 12.587587;  WData[4].Value[6] := 366;   WData[4].Value[7] := 5;   WData[4].Value[8] := 16;  WData[4].Value[9] := 11;
      WData[4].Value[10] := 4.94874;  WData[4].Value[11] := 4445.777;
      WData[0].Name := 'Ht'; WData[1].Name := 'abracadabra'; WData[2].Name := 'toolonganame'; WData[3].Name := 'A089478464764'; WData[4].Name := 'baxter56';
      writeln('Data created successfully. 12 x 5. Big data. Vars 4 and 5 correlated. Two data are NaN.');
    end;
    'T2': begin
      nRow := 2;
      nCol := 2;
      SetLength(WData, nCol);   //invertible
      for i := 0 to nCol - 1 do
        SetLength(WData[i].Value, nRow);
      WData[0].Name := 'First'; WData[1].Name := 'Second';
      WData[0].Value[0] := 1; WData[0].Value[1] := 0;
      WData[1].Value[0] := 0; WData[1].Value[1] := 2;
      writeln('Data created successfully. 2 x 2.');
    end;
    'T3': begin
      nRow := 3;
      nCol := 3;
      SetLength(WData, nCol);             //invertible det =-9
      for i := 0 to nCol - 1 do
        SetLength(WData[i].Value, nRow);
      WData[0].Name := 'Way'; WData[1].Name := 'Here'; WData[2].Name := 'Avenue';
      WData[0].Value[0] := 2; WData[0].Value[1] := 1; WData[0].Value[2] := 3;
      WData[1].Value[0] := 0; WData[1].Value[1] := 1; WData[1].Value[2] := 4;
      WData[2].Value[0] := 5; WData[2].Value[1] := 2; WData[2].Value[2] := 1;
      writeln('Data created successfully. 3 x 3 invertible.');
    end;
    'T4': begin
     nRow := 5;
     nCol := 3;
     SetLength(WData, nCol);
     for i := 0 to nCol - 1 do
       SetLength(WData[i].Value, nRow);
      WData[0].Name := 'St'; WData[1].Name := 'Yeah'; WData[2].Name := 'Avenue';
      WData[0].Value[0] := 2.1;   WData[0].Value[1] := 2.2;   WData[0].Value[2] := 2.3;   WData[0].Value[3] := 2.4;   WData[0].Value[4] := 2.5;
      WData[1].Value[0] := 3.4;   WData[1].Value[1] := 3.51;   WData[1].Value[2] := 3.62;   WData[1].Value[3] := 3.77;   WData[1].Value[4] := 3.86;
      WData[2].Value[0] := 3;   WData[2].Value[1] := 4;   WData[2].Value[2] := 9;   WData[2].Value[3] := 2;   WData[2].Value[4] := 1;
      writeln('Data created successfully. 5 x 3.');
    end;
    'T5': begin
      nRow := 2;
      nCol := 2;
      SetLength(WData, nCol);   //invertible
      for i := 0 to nCol - 1 do
        SetLength(WData[i].Value, nRow);
      WData[0].Name := 'First000000'; WData[1].Name := 'Second00000000';
      WData[0].Value[0] := 0; WData[0].Value[1] := 0;
      WData[1].Value[0] := 0; WData[1].Value[1] := 0;
      writeln('Data created successfully. 2 x 2 of zeroes.');
    end;
    'T6': begin
      nRow := 2;
      nCol := 2;
      SetLength(WData, nCol);   //invertible
      for i := 0 to nCol - 1 do
        SetLength(WData[i].Value, nRow);
      WData[0].Name := 'First0'; WData[1].Name := 'Second0';
      WData[0].Value[0] := 1; WData[0].Value[1] := 1;
      WData[1].Value[0] := 1; WData[1].Value[1] := 1;
      writeln('Data created successfully. 2 x 2 of ones.');
    end;
    'T7': begin
      nRow := 2;
      nCol := 1;
      SetLength(WData, nCol);   //invertible
      for i := 0 to nCol - 1 do
        SetLength(WData[i].Value, nRow);
      WData[0].Name := 'First0';
      WData[0].Value[0] := 2; WData[0].Value[1] := 4;
      writeln('Data created successfully. 2 x 1 of [2, 4].');
    end;
    'T8': begin
      nRow := 1;
      nCol := 1;
      SetLength(WData, nCol);   //invertible
      for i := 0 to nCol - 1 do
        SetLength(WData[i].Value, nRow);
      WData[0].Name := 'Firstofall';
      WData[0].Value[0] := 0;
      // Process NaN.
      ConvertNaNinData;
      //Fill the caches for WData.
      FillWDataCaches;
      writeln('Data created successfully. 1 x 1 of zero.');
    end;
    'T10': begin
      nRow := 4;
      nCol := 4;
      SetLength(WData, nCol);   // Invertible, determinant = 87.
      for i := 0 to nCol - 1 do
        SetLength(WData[i].Value, nRow);
      // Column names
      WData[0].Name := 'Way';
      WData[1].Name := 'Here';
      WData[2].Name := 'Avenue';
      WData[3].Name := 'Road';

      // Column 0 (Way)
      WData[0].Value[0] := 4;
      WData[0].Value[1] := 0;
      WData[0].Value[2] := 7;
      WData[0].Value[3] := 1;

      // Column 1 (Here)
      WData[1].Value[0] := 1;
      WData[1].Value[1] := 5;
      WData[1].Value[2] := 0;
      WData[1].Value[3] := 2;

      // Column 2 (Avenue)
      WData[2].Value[0] := 3;
      WData[2].Value[1] := 2;
      WData[2].Value[2] := 6;
      WData[2].Value[3] := 0;

      // Column 3 (Road)
      WData[3].Value[0] := 2;
      WData[3].Value[1] := 1;
      WData[3].Value[2] := 4;
      WData[3].Value[3] := 3;
      writeln('Data created successfully. 4 x 4 invertible matrix with determinant of 87.');
    end;
    'T11': begin
      nRow := 3;
      nCol := 3;
      SetLength(WData, nCol);             //invertible det =-9
      for i := 0 to nCol - 1 do
        SetLength(WData[i].Value, nRow);
      WData[0].Name := 'Way'; WData[1].Name := 'Here'; WData[2].Name := 'Avenue';
      WData[0].Value[0] := 2; WData[0].Value[1] := 1; WData[0].Value[2] := 0;
      WData[1].Value[0] := 1; WData[1].Value[1] := 2; WData[1].Value[2] := 1;
      WData[2].Value[0] := 0; WData[2].Value[1] := 1; WData[2].Value[2] := 2;
      writeln('Data created successfully. 3 x 3 invertible.');
    end;
  end;
  // Process NaN.
  ConvertNaNinData;
  // Fill the caches for WData.
  FillWDataCaches;
end;

// For a new or modified variable, fill the caches.
procedure FillVectorCaches(var Col: WVectorType);
begin
  // Variable without name gets "Var" plus a number.
  if Col.Name = '' then begin
    Col.Name := 'Var' + IntToStr(TVN);
    Inc(TVN);
  end;
  Col.Median := Median(Col);
  Col.Mean := ArithmeticMean(Col);
  Col.StdDev := StandardDeviation(Col, True);
  Col.Min := MinmValue(Col);
  Col.Max := MaxmValue(Col);
end;

// Fill caches of WData all at once. Used in test data and manually-input data.
procedure FillWDataCaches;
Var
  i: Integer;
begin
  for i := 0 to High(WData) do
    FillVectorCaches(WData[i]);
  SuppressErrorStack;  //suppress error messages
end;

{ Access Data }

// Put data into WData directly.
procedure PutWVector(const Col: WVectorType; const iCol: Integer);
begin
  DeepCopyVector(Col, WData[iCol]);
end;

// Get data from WData directly.
procedure GetWVector(var Col: WVectorType; const iCol: Integer);
begin
  DeepCopyVector(WData[iCol], Col);
end;

// Add a new (or replace a) vector.
procedure AddWVector(const NewName: string; const NewVar: WVectorType; var AddError: Boolean);
var
  i, m: Integer;
begin
  AddError := False;
  m := Length(NewVar.Value);

  // Check if already variable name, and replace it.
  for i := 0 to nCol - 1 do  begin
    if UpCase(WData[i].Name) = UpCase(NewName) then begin
      DeepCopyVector(NewVar, WData[i]);
      FillVectorCaches(WData[i]);
      Exit;
    end;
  end;

  // Check enough room.
{  if nCol >= 99 then begin
    writeln('Not enough space for variable.');
    AddError := True;
    Exit;
  end; Add try except}

  // If not found above, then add a new variable.
  // Proper name has already been checked, so add variable, and fill.
  SetLength(WData, nCol + 1);
  SetLength(WData[nCol].Value, m);
  SetLength(WData[nCol].IsNaN, m);
  // Copy variable.
  DeepCopyVector(NewVar, WData[nCol]);
  WData[nCol].Name := NewName;
  // Set IsNan to False.
{  writeln('2add var nrow m ncol', nrow, ' ', m, ' ', ncol);
  for j := 0 to m - 1 do begin
    writeln(j);
    WData[nCol].IsNaN[j] := False;
  end;}
  // Fill the statistics caches.         Alreadt filled"
  FillVectorCaches(WData[nCol]);
  Inc(nCol);
end;

// Add a new (or replace a) scalar/constant.
procedure AddWConstant(const NewName: string; const NewVar: WVectorType; var AddError: Boolean);
var
  i: Integer;
begin
  AddError := False;
  for i := 0 to nScal - 1 do begin
    // Already constant with this name.
    if UpCase(SData[i].Name) = UpCase(NewName) then begin
      SData[i].Value := NewVar.Value[0];
      Exit;
    end;
  end;

  // If not found , then add a new new scalar.
  // Proper name has already been checked, so add constant.
  SetLength(SData, nScal + 1);
  SData[nScal].Value := NewVar.Value[0];
  SData[nScal].Name := NewName;
  Inc(nScal);
end;

// Fill SData with some nice, precalculated scalars.
procedure FillScalars;
begin
  SetLength(SData, nScal + 15);
  SData[nScal].Name :=  'Pi'; SData[nScal].Value :=  Pi;
  SData[nScal + 1].Name := 'TwoPi'; SData[nScal + 1].Value := TwoPi;
  SData[nScal + 2].Name := 'InvSqrt2'; SData[nScal + 2].Value := InvSqrt2;
  SData[nScal + 3].Name :=  'SqrtPi2'; SData[nScal + 3].Value := SqrtPi2;
  SData[nScal + 4].Name :=  'Ln2'; SData[nScal + 4].Value := Ln2;
  SData[nScal + 5].Name :=  'Sqrt2Pi'; SData[nScal + 5].Value := Sqrt2Pi;
  SData[nScal + 6].Name :=  'Sqrt2'; SData[nScal + 6].Value := Sqrt2;
  SData[nScal + 7].Name :=  'Sqrt3'; SData[nScal + 7].Value := Sqrt2;
  SData[nScal + 8].Name :=  'SD1'; SData[nScal + 8].Value := SD1;
  SData[nScal + 9].Name :=  'SD2'; SData[nScal + 9].Value := SD2;
  SData[nScal + 10].Name :=  'SD3'; SData[nScal + 10].Value := SD3;
  SData[nScal + 11].Name :=  'E'; SData[nScal + 11].Value := E;
  SData[nScal + 12].Name :=  'Ln10'; SData[nScal + 12].Value := Ln10;
  SData[nScal + 13].Name :=  'Delta'; SData[nScal + 13].Value := Delta;
  SData[nScal + 14].Name :=  'Phi'; SData[nScal + 14].Value := Phi;
  Inc(nScal, 14);
end;

{ Read & Save Data }
procedure ReadCSVData(const FileName: String; var RowCount, ColCount: Integer);
var
  F: TextFile;
  Line: String;
  Tokens: SVectorType;
  TempTable: CMatrixType;  // temporary: obs × vars
  FirstRow: Boolean;
  i, j, ObsIndex, Code: Integer;
  Handle: THandle;
begin
  RowCount := 0;    // number of observations
  ColCount := 0;    // number of variables
  FirstRow := True;

  // Check file availability.
  if not FileExists(FileName) then begin
    writeln('File not found: ', FileName);
    Exit;
  end;
  Handle := FileOpen(FileName, fmOpenRead or fmShareExclusive);
  if Handle = -1 then begin
    writeln('File is already open elsewhere.');
    Exit;
  end;

  // Try opening file.
  try
    AssignFile(F, FileName);
    Reset(F);

    // Read CSV into Temp structure.
    ObsIndex := 0;
    while not EOF(F) do begin
      ReadLn(F, Line);
      Tokens := Line.Split([',']);

      // First line consists of variable names.
      if FirstRow then begin
        FirstRow := False;
        ColCount := Length(Tokens);

        // Store names.
        SetLength(WData, ColCount);
        for i := 0 to ColCount - 1 do
          WData[i].Name := Tokens[i];

        // Initialize temp table with zero observations.
        SetLength(TempTable, 0);
        continue;
      end;

      // Check column count.
      if Length(Tokens) <> ColCount then begin
        writeln('Column count mismatch at CSV row ', ObsIndex + 2);
        CloseFile(F);
        Exit;
      end;

      // Extend temp table (add observation).
      Inc(RowCount);
      SetLength(TempTable, RowCount);
      SetLength(TempTable[RowCount - 1], ColCount);

      // Convert strings → numbers.
      for j := 0 to ColCount - 1 do begin
        Val(Tokens[j], TempTable[RowCount - 1][j], Code);
        if Code <> 0 then begin
          if Verbose or DebugOn then
            writeln('Conversion error at row ', RowCount + 1,
              ' col ', j + 1, '; value "', Tokens[j], '"; input as NaN.');
          TempTable[RowCount - 1, j] := NaN;
        end;
      end;
    end;
    CloseFile(F);

  except
    writeln('File cannot be accessed.');
  end;

  // Transpose data to WData format.
  for j := 0 to ColCount - 1 do begin
    // RowCount observations per variable.
    SetLength(WData[j].Value, RowCount);
    SetLength(WData[j].IsNan, RowCount);
  end;

  // Deal with NaN.
  ConvertNaNinData;

  // Fill the statistics caches.
  FillWDataCaches;

  // Add flag if each value was a NaN value may be changed later).
  for i := 0 to RowCount - 1 do
    for j := 0 to ColCount - 1 do begin
      WData[j].Value[i] := TempTable[i, j];
      if IsNan(TempTable[i, j]) then
        WData[j].IsNan[i] := True
      else
        WData[j].IsNan[i] := False;
    end;

  writeln('Data read successfully from ', FileName, '. There are ',
    nCol, ' variables, and ', nRow, ' observations.');
end;

procedure SavePartialData(const Mess: String; PProc: ProcType);  overload;
var
  WInput: Char;
  FileName: String;
  SavedOutput: Text;
begin
  FirstPass := False;
  write(Mess);
  readln(WInput);
  if UpCase(WInput) = 'Y' then begin
    write('Name of file: ');
    read(FileName);

    // Save console output
    SavedOutput := Output;

    // Redirect to file.
    Flush(Output);
    Assign(Output, FileName);
    Rewrite(Output);
    PProc;

    // Restore console.
    Flush(Output);
    writeln;
    Output := SavedOutput;
    writeln('Save complete to ', FileName, '.');
  end;
end;

procedure SavePartialData(const Mess: String; const Col: WVectorType; VProc: VProcType);  overload;
var
  WInput: Char;
  FileName: String;
  SavedOutput: Text;
begin
  FirstPass := False;
  write(Mess);
  readln(WInput);
  if UpCase(WInput) = 'Y' then begin
    write('Name of file: ');
    read(FileName);

    // Save console output
    SavedOutput := Output;

    // Redirect to file.
    Flush(Output);
    Assign(Output, FileName);
    Rewrite(Output);
    VProc(Col);

    // Restore console.
    Flush(Output);
    writeln;
    Output := SavedOutput;
    writeln('Save complete to ', FileName, '.');
  end;
end;

procedure SavePartialData(const Mess: String; const Col1, Col2: WVectorType; VVProc: VVProcType);  overload;
var
  WInput: Char;
  FileName: String;
  SavedOutput: Text;
begin
  FirstPass := False;
  write(Mess);
  readln(WInput);
  if UpCase(WInput) = 'Y' then begin
    write('Name of file: ');
    read(FileName);

    // Save console output
    SavedOutput := Output;

    // Redirect to file.
    Flush(Output);
    Assign(Output, FileName);
    Rewrite(Output);
    VVProc(Col1, Col2);

    // Restore console.
    Flush(Output);
    writeln;
    Output := SavedOutput;
    writeln('Save complete to ', FileName, '.');
  end;
end;

procedure SavePartialData(const Mess: String; const WVar: IVectorType; IProc: IProcType);  overload;
var
  WInput: Char;
  FileName: String;
  SavedOutput: Text;
begin
  FirstPass := False;
  write(Mess);
  readln(WInput);
  if UpCase(WInput) = 'Y' then begin
    write('Name of file: ');
    read(FileName);

    // Save console output
    SavedOutput := Output;

    // Redirect to file.
    Flush(Output);
    Assign(Output, FileName);
    Rewrite(Output);
    IProc(WVar);

    // Restore console.
    Flush(Output);
    writeln;
    Output := SavedOutput;
    writeln('Save complete to ', FileName, '.');
  end;
end;

procedure SaveData(const FileName: String);
var
  F: TextFile;
  i, j: Integer;
  RowStr: String;
  NumVars, NumObs: Integer;
begin
  NumVars := Length(WData);
  if NumVars = 0 then begin
    writeln('No data to save.');
    Exit;
  end;

  NumObs := Length(WData[0].Value);
  if NumObs = 0 then begin
    writeln('No data to save.');
    Exit;
  end;

  AssignFile(F, FileName);
  Rewrite(F);

  // Write header.
  RowStr := '';
  for j := 0 to NumVars - 1 do begin
    if j > 0 then RowStr := RowStr + ',';
    RowStr := RowStr + WData[j].Name;
  end;
  Writeln(F, RowStr);

  // Write observations.
  for i := 0 to NumObs - 1 do begin
    RowStr := '';
    for j := 0 to NumVars - 1 do begin
      if j > 0 then RowStr := RowStr + ',';
      RowStr := RowStr + FloatToStr(WData[j].Value[i]);
    end;
    Writeln(F, RowStr);
  end;

  CloseFile(F);
  writeln('Data saved successfully as ', FileName, '.');
end;

procedure SaveFile(const X: CMatrixType; const FileName: String); overload;
var
  F: TextFile;
  i, j: Integer;
  RowStr: String;
  NumVars, NumObs: Integer;
begin
  NumVars := Length(X);
  if NumVars = 0 then begin
    writeln('No data to save.');
    Exit;
  end;

  NumObs := Length(X[0]);
  if NumObs = 0 then begin
    writeln('No data to save.');
    Exit;
  end;

  AssignFile(F, FileName);
  Rewrite(F);

  // Write observations.
  for i := 0 to NumObs - 1 do begin
    RowStr := '';
    for j := 0 to NumVars - 1 do begin
      if j > 0 then RowStr := RowStr + ',';
      RowStr := RowStr + FloatToStr(X[j, i]);
    end;
    Writeln(F, RowStr);
  end;

  CloseFile(F);
  writeln('Data saved successfully as ', FileName, '.');
end;

procedure SaveFile(const X: WMatrixType; const FileName: String); overload;
var
  F: TextFile;
  i, j: Integer;
  RowStr: String;
  NumVars, NumObs: Integer;
begin
  NumVars := Length(X);
  if NumVars = 0 then begin
    writeln('No data to save.');
    Exit;
  end;

  NumObs := Length(X[0].Value);
  if NumObs = 0 then begin
    writeln('No data to save.');
    Exit;
  end;

  AssignFile(F, FileName);
  Rewrite(F);

  // Write observations.
  for i := 0 to NumObs - 1 do begin
    RowStr := '';
    for j := 0 to NumVars - 1 do begin
      if j > 0 then RowStr := RowStr + ',';
      RowStr := RowStr + FloatToStr(X[j].Value[i]);
    end;
    Writeln(F, RowStr);
  end;

  CloseFile(F);
  writeln('Data saved successfully as ', FileName, '.');
end;

function TeeWrite(var F: TextRec): Integer;
var
  S: String;
begin
  // Extract buffered text.
  SetString(S, PChar(F.BufPtr), F.BufPos);

  // Write directly to console.
  if S <> '' then
    FileWrite(TextRec(Output).Handle, S[1], Length(S));

  // Write to tee file.
  if TeeActive and (S <> '') then
    FileWrite(TextRec(TeeFile).Handle, S[1], Length(S));

  F.BufPos := 0;
  Result := 0;
end;

procedure StartTeeOutput(const FileName: String);
begin
  writeln('Start logging to ', FileName, '.');
  TeeFileName := FileName;
  AssignFile(TeeFile, FileName);

  {$I-}
  Append(TeeFile);          // Try to append.
  if IOResult <> 0 then
    Rewrite(TeeFile);       // If append failed, create new file.
  {$I+}

  TeeActive := True;

  // Save original handlers.
  OldInOut  := TextRec(Output).InOutFunc;
  OldFlush  := TextRec(Output).FlushFunc;

  // Install tee hook.
  with TextRec(Output) do begin
    Mode      := fmOutput;
    BufPos    := 0;
    InOutFunc := nil;
    FlushFunc := @TeeWrite;
  end;
end;

{procedure StartTeeOutput(const FileName: String);
begin
  TeeFileName := FileName;

  AssignFile(TeeFile, FileName);
  Rewrite(TeeFile);
  TeeActive := True;

  // Save original handlers.
  OldInOut  := TextRec(Output).InOutFunc;
  OldFlush  := TextRec(Output).FlushFunc;

  // Install tee hook.
  with TextRec(Output) do
  begin
    Mode     := fmOutput;
    BufPos   := 0;
    InOutFunc := nil;             // Disable normal write.
    FlushFunc := @TeeWrite;       // Custom flush.
  end;
end;}

procedure StopTeeOutput;
begin
  if TeeActive then
  begin
    TeeActive := False;

    // Restore original handlers.
    TextRec(Output).InOutFunc := OldInOut;
    TextRec(Output).FlushFunc := OldFlush;

    Flush(Output);                // Flush console normally.
    CloseFile(TeeFile);           // Close log file.

    Writeln('Logging stopped to ', TeeFileName);
  end;
end;

{  Modify Routines }
procedure DeleteVariable(const iCol: Integer);  // Zero based.
var
  i: Integer;
  TempName: String;
begin
  TempName := WData[iCol].Name;
  for i := iCol to nCol - 2 do
    DeepCopyVector(WData[i + 1], WData[i]);
  Dec(nCol);
  SetLength(WData, Length(WData) - 1);
  writeln('Variable ', iCol + 1, ' ', TempName, ' has been deleted.');
end;

procedure DeleteScalar(const iScal: Integer);    // Zero based.
var
  i: Integer;
  TempName: String;
begin
  TempName := SData[iScal].Name;
  for i := iScal to nScal - 2 do
    SData[i] := SData[i + 1];
  Dec(nScal);
  SetLength(SData, Length(SData) - 1);
  writeln('Constant ', iScal + 1, ' ', TempName, ' has been deleted.');
end;

procedure DeleteRow(const iRow: Integer);        // Zero based.
var
  i, j: Integer;
begin
  for i := 0 to nCol - 1 do begin
    for j := iRow to nRow - 2 do
      WData[i].Value[j] := WData[i].Value[j + 1];
    SetLength(WData[i].Value, Length(WData[i].Value) - 1);
  end;
  Dec(nRow);
  writeln('Row ', iRow + 1, ' has been deleted.');
end;

procedure DataTrim(const iCol: Integer; const L, R: FloatType);
var
  i, j, k, n, m: Integer;
begin                        // Use another var for n.
  n := Length(WData);
  for j := 0 to nRow - 1 do
    if (WData[iCol].Value[j] <= L) or (WData[iCol].Value[j] >= R) then
      for i := 0 to n - 1 do begin
        m := Length(WData[i].Value);
        for k := i to  m - 1 do
          WData[i].Value[k] := WData[i].Value[k + 1];
        SetLength(WData[i].Value, m - 1);
      end;
  for i:= 1 to n - 1 do
    if Length(WData[i].Value) > nRow then nRow := Length(WData[i].Value);
  if nRow <= 0 then
    writeln('Error on trimming: all rows trimmed.');
end;

procedure WinsorizeWColumn(var Col: WVectorType; LowerPct, UpperPct: FloatType);
var
  LowerBound, UpperBound: FloatType;
  j, m: Integer;
  SortedCol: WVectorType;
begin
  DeepCopyVector(Col, SortedCol);
  m := Length(SortedCol.Value);

  // Sort.
  QuickSort(SortedCol, 0, m);

  LowerBound := SortedCol.Value[Round(LowerPct * m)];
  UpperBound := SortedCol.Value[Round(UpperPct * m)];

  // Winsorize original data.
  for j := 0 to m - 1 do begin
    if SortedCol.Value[j] < lowerBound then
      SortedCol.Value[j] := LowerBound
    else if Col.Value[j] > upperBound then
      SortedCol.Value[j] := upperBound;
  end;

  //Return vector to WData.
  FillVectorCaches(SortedCol);
  DeepCopyVector(SortedCol, Col);
end;

procedure WinsorizeColumn(var Col: WVectorType; iCol: Integer; LowerPct, UpperPct: FloatType);
var
  LowerBound, UpperBound: FloatType;
  j, m: Integer;
  SortedCol: WVectorType;
begin
  // Get the vector
  DeepCopyVector(Col, SortedCol);
  m := Length(SortedCol.Value);

  // Sort
  QuickSort(SortedCol, 0, m);

  LowerBound := SortedCol.Value[Round(LowerPct * m)];
  UpperBound := SortedCol.Value[Round(UpperPct * m)];

  // Winsorize original data
  for j := 0 to m - 1 do begin
    if SortedCol.Value[j] < lowerBound then
      SortedCol.Value[j] := LowerBound
    else if SortedCol.Value[j] > upperBound then
      SortedCol.Value[j] := upperBound;
  end;

  // Return vector to WData
  FillVectorCaches(SortedCol);
  PutWVector(SortedCol, iCol);
end;

procedure YeoJohnsonTransform(var Data: WVectorType; Lambda: FloatType);
var
  j, m: Integer;         // Rename data as ccol.
  xi: FloatType;
begin
  m := Length(Data.Value);
  for j := 0 to m - 1 do begin
    xi := Data.Value[j];
    if xi >= 0 then begin
      if Lambda = 0 then
        Data.Value[j] := Ln(xi + 1)
      else
        Data.Value[j] := (Power(xi + 1, Lambda) - 1) / Lambda;
    end
    else begin
      if Lambda = 2 then
        Data.Value[j] := -Ln(-xi + 1)
      else
        Data.Value[j] := - (Power(-xi + 1, 2 - Lambda) - 1) / (2 - Lambda);
    end;
  end;
  FillvectorCaches(Data);
end;

{Centering and Standardizing Routines}
// Center the column by making mean = 0.
procedure CenterColumn(var Col: WVectorType); overload;
var
  mu: FloatType;
  j, m: Integer;
begin
  m := Length(Col.Value);
  mu := ArithmeticMean(Col);
  for j := 0 to m - 1 do
    Col.Value[j] := Col.Value[j] - mu;
  FillVectorCaches(Col);
end;

procedure CenterColumn(var Col: CVectorType); overload;
var
  mu: FloatType;
  j, m: Integer;
begin
  m := Length(Col);
  mu := Mean(Col);
  for j := 0 to m - 1 do
    Col[j] := Col[j] - mu;
end;

//Center the column by making mean = 0, and make SD = 1.
procedure StandardizeColumn(var Col: WVectorType); overload;
var
  mu, Sigma: FloatType;
  j, m: Integer;
begin
  m := Length(Col.Value);
  mu := ArithmeticMean(Col);
  sigma := Sqrt(SumOfSquares(Col) / (m - 1));
  for j := 0 to m - 1 do begin
    if (sigma <> 0.0) then
      Col.Value[j] := (Col.Value[j] - mu) / sigma
    else
      Col.Value[j] := 0.0;
  end;
  FillVectorCaches(Col);
end;

procedure UnstandardizeColumn(var Col: WVectorType; const mu, Sigma: FloatType);
var
  j, m: Integer;
begin
  m := Length(Col.Value);
  for j := 0 to m - 1 do
    Col.Value[j] := (Col.Value[j] * sigma) + mu;
end;

procedure StandardizeRetainWColumn(var Col: WVectorType; var mu, Sigma: FloatType);
var
  j, m: Integer;
begin
  m := Length(Col.Value);
  mu := ArithmeticMean(Col);
  sigma := Sqrt(SumOfSquares(Col) / (m - 1));
  for j := 0 to m - 1 do begin
    if (sigma <> 0.0) then
      Col.Value[j] := (Col.Value[j] - mu) / sigma
    else
      Col.Value[j] := 0.0;
  end;
end;

procedure StandardizeColumn(var Col: CVectorType); overload;
var
  mu, Sigma: FloatType;
  j, m: Integer;
begin
  m := Length(Col);
  mu := Mean(Col);
  sigma := Sqrt(SumOfSquares(Col) / (m - 1));
  for j := 0 to m - 1 do begin
    if (sigma <> 0.0) then
      Col[j] := (Col[j] - mu) / sigma
    else
      Col[j] := 0.0;
  end;
end;

procedure CenterMatrix(const Data: WMatrixType);
var
  i, n: Integer;
begin
  n := Length(Data);
  for i := 0 to n - 1 do
    CenterColumn(Data[i]);
end;

procedure StandardizeMatrix(const Data: WMatrixType); overload;
var
  i, n: Integer;
begin
  n := Length(Data);
  for i := 0 to n - 1 do
    StandardizeColumn(Data[i]);
end;

procedure StandardizeMatrix(const Data: CMatrixType); overload;
var
  i, n: Integer;
begin
  n := Length(Data);
  for i := 0 to n - 1 do
    StandardizeColumn(Data[i]);
end;

end.
