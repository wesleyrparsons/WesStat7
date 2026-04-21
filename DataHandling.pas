unit DataManager;

{$mode ObjFPC}{$H+}

{ WesStat7 September 24, 2025
  Source code by Wesley R. Parsons,
  wespar@bellouth.net, www.wespar.com, Miami, Florida, 2025.}

interface

uses
  BasicFuncs,
  Crt,
  Globals,
  Math,
  SysUtils;

function SmartFloat(const Value: Extended): String;
function WriteTab(const Txt: String; const Value: Extended): String;
function WriteTabln(const Txt: String; const Value: Extended): String;
function WriteTabPeriodln(const Txt: String; const Value: Extended): String;
function CardinalToOrdinal(const n: Integer): String;
procedure CheckWrite(Col: TVector);
function ReadPositiveInteger(const Pr: string; maxl: Integer): integer;
procedure CreateRandomData(var Data: TMatrix);
procedure DataCorrelateAndZero(var Data: TMatrix);
procedure DisplayCorrelationMatrix(const CData, SData: TCMatrix; const WVar: TIVector;
  const IsPearson, PrintSig: Boolean; const CorrMess: TCorr);
procedure DisplayCovarianceMatrix(const CVData: TCMatrix; const WVar: TIVector);
procedure DisplaySelectedData(const Data: TMatrix; const Condensed: Boolean);
procedure DisplayData(const PData: TMatrix; const PLabel: TIVector);
procedure DisplayDataA(const PData: TMatrix; WVar: TString);
procedure DisplayDataCondensed(const PData: TMatrix; const PLabel: TIVector);
procedure DisplayMatrix(const PData: TMatrix);
procedure DisplayMatrixM(const Mess: String; const PData: TMatrix);
procedure DisplayMatrixAll;
procedure DC(const WVar: TIVector);
procedure DME(const Mess: String; const Data: TMatrix);
procedure DE(const Data: TMatrix);
procedure DisplayColumn(const PData: TVector);
procedure DisplayColumnM(const Mess: String; const PData: TVector);
procedure SaveData(const Mess: String; const Data: TMatrix);
procedure ReadCSVData(const FileName: string; var M: TMatrix; var RowCount, ColCount: Integer);
procedure SelectColumnsA(const Data: TMatrix; var VData: TMatrix;
  var WVar: TIVector; var nVar: Integer);
procedure ManuallyInputData(var WData: TMatrix; var nRow, nCol: Integer);
procedure ModifyData(var WData: TMatrix);
procedure GetTestData(const Str: String);

implementation

function SmartFloat(const Value: Extended): String;
var
  //Pad: Integer;
  Num: String;
begin
  if IsNan(Value) then
    Num := 'NaN'
  else if (Value = 0) then Num := '0'
  else if (abs(Value) < 1e-4) or (abs(Value) > 1e6) then
    Num := FormatFloat('0.000E+00', Value)  // Scientific notation
  else
    Num := FormatFloat('#,##0.####', Value); // Regular format with up to 4 decimals
  //Pad := Length(Num);
  Result := Num;
  //Result := Num + StringofChar(' ', 40 - Pad);   //no padding, use WriteTab for that
end;

function WriteTab(const Txt: String; const Value: Extended): String;
var
  Pad: Integer;
  Num: String;
begin
  if IsNan(Value) then
    Num := 'NaN'
  else if (Value = 0) then
    Num := '0'
   // else
  else if (abs(Value) < 1e-4) or (abs(Value) > 1e6) then
    Num := FormatFloat('0.0000E+00', Value)  // Scientific notation
  else
    Num := FormatFloat('#,##0.#####', Value); // Regular format with up to 5 decimals
  Pad := Length(Num) + Length(Txt);
  if Pad > 44 then Pad := 45
  else
    Result := Txt + Num + StringofChar(' ', 45 - Pad);
  Write(Result);
end;

function WriteTabln(const Txt: String; const Value: Extended): String;
var
  Num: String;
begin
   if IsNan(Value) then
    Num := 'NaN'
  else if (Value = 0) then Num := '0'
  else if (abs(Value) < 1e-4) or (abs(Value) > 1e6) then
    Num := FormatFloat('0.0000E+00', Value)  // Scientific notation
  else
    Num := FormatFloat('#,##0.#####', Value); // Regular format with up to 5 decimals
  Result := Txt + Num;
  writeln(Result);
end;

function WriteTabPeriodln(const Txt: String; const Value: Extended): String;
var
  Num: String;
begin
   if IsNan(Value) then
    Num := 'NaN'
  else if (Value = 0) then Num := '0'
  else if (abs(Value) < 1e-4) or (abs(Value) > 1e6) then
    Num := FormatFloat('0.000E+00', Value)  // Scientific notation
  else
    Num := FormatFloat('#,##0.####', Value); // Regular format with up to 4 decimals
  Result := Txt + Num;
  writeln(Result, '.');
end;

function CardinalToOrdinal(const n: Integer): String;
var
  suffix: String;
  lastTwo, lastDigit: Integer;
begin
  lastTwo := n mod 100;
  lastDigit := n mod 10;
  // Handle special cases: 11th, 12th, 13th
  if (lastTwo >= 11) and (lastTwo <= 13) then
    suffix := 'th'
  else
    case lastDigit of
      1: suffix := 'st';
      2: suffix := 'nd';
      3: suffix := 'rd';
    else
      suffix := 'th';
    end;
  Result := IntToStr(n) + suffix;
end;

procedure CheckWrite(Col: TVector);
var
  j: Integer;
begin
  writeln('Check Column');
  for j := 0 to High(Col.Value) do Write(j, '  ', Col.Value[j], '     ');
  writeln;
end;

function ReadPositiveInteger(const Pr: string; maxl: integer): integer;
var
  InputStr: string;
  Value, Code: integer;
begin
  repeat
    Write(Pr);
    ReadLn(InputStr);
    Val(InputStr, Value, Code);  // try to convert string to integer
    if (Code <> 0) or (Value <= 0) then
      Writeln('Please enter a valid whole number greater than 0.');
    if Value > maxl then
    begin
      Writeln('Please enter a number equal to or less than ', maxl, '.');
      Code := -99;
    end;
  until (Code = 0) and (Value > 0);
  Result := Value;
end;

function RandomNormal(const mu, sigma: TFloat): TFloat;
var
  u1, u2: TFloat;
begin
  repeat
    u1 := Random;
  until u1 > 0.0;
  u2 := Random;
  Result := mu + sigma * sqrt(-2.0 * ln(u1)) * cos(2.0 * PI * u2);
end;

procedure CreateRandomData(var Data: TMatrix);
var
  i, j: word;
begin
  for i := 0 to High(Data) do
    for j := 0 to High(Data[0].Value) do
      Data[i].Value[j] := RandomNormal(10, 5);
end;

procedure DataCorrelateAndZero(var Data: TMatrix);
var
  i, j: Integer;
begin
  for i := 0 to High(Data) do
    for j := 0 to High(Data[0].Value) do begin
      Data[2].Value[j] := abs(Data[2].Value[j]) + 5.0;  //  col 3 > 0
      Data[2].Value[j] := Data[0].Value[j] + Random;  //highly correlated
    end;
end;

procedure ManuallyInputData(var Data: TMatrix);
var
  WDataInput: TFloat;
  i, j: Word;
begin
  Writeln('Enter each piece of Data, follow by <CR>.');
  for j := 0 to High(Data) do
  begin
    Write('Row ', j);
    for i := 0 to High(Data[j].Value) do
    begin
      Read(WDataInput);
      Data[i].Value[j] := WDataInput;
    end;
    writeln;
  end;
end;

procedure ReadCSVData(const FileName: string; var M: TMatrix; var RowCount, ColCount: Integer);
var
  F: TextFile;
  Line: string;
  Tokens: TStringArray;
  TempRow: array of TFloat;
  i, j: Integer;
begin
  RowCount := 0;
  ColCount := 0;
  if FileExists(FileName) then begin
    AssignFile(f, FileName);
    Reset(f);
  end
  else
    writeLn('File not found: ', FileName);
  SetLength(M, 0); // Start with empty matrix
  while not EOF(F) do begin
    ReadLn(F, Line);
    Tokens := Line.Split([',']);
    if ColCount = 0 then
      ColCount := Length(Tokens)
    else if Length(Tokens) <> ColCount then begin
      writeln('Error: Inconsistent column count at row ', RowCount + 1);
      writeln('File not read.');
      Exit;
    end;
    SetLength(TempRow, ColCount);
    for i := 0 to ColCount - 1 do
      Val(Tokens[i], TempRow[i], j); // j is error code; 0 = success
      if j <> 0 then begin
        writeln('Error: Data not read on row ', i + 1, ', column ', j + 1);
      end;
    Inc(RowCount);
    SetLength(M, RowCount);
    M[RowCount - 1].Value := TempRow;
  end;
  CloseFile(F);
end;

procedure DisplayCorrelationMatrix(const CData, SData: TCMatrix; const WVar: TIVector;
  const IsPearson, PrintSig: Boolean; const CorrMess: TCorr);
var
  i, j, n: integer;
begin
  n := Length(CData);
  Case CorrMess of
    Pearson: write('Pearson''s R');
    Spearman: write('Spearman''s Rho');
    Kendall: write('Kendall''s Tau');
    Hoeffding: write('Hoeffding'' D');
  end;
  writeln(' Correlation Matrix   Variables: ', n);
  write('            ');
  for j := 0 to n - 1 do
    Write('Var', WVar[j]: 3, '    ');
  writeln;
  for i := 0 to n - 1 do begin
    Write('Var', WVar[i]: 3, '    ');
    for j := 0 to n - 1 do Write(CData[i, j]: 10: 6);
    writeln;
    if PrintSig then begin
      Write('          ');
      for j := 0 to n - 1 do
        Write(' (', SData[i, j]: 6: 5, ')');
      writeln;
    end;
  end;
  if PrintSig then writeln('Two-tailed significance in parentheses.');
  Pause;
  writeln;
end;

procedure DisplayCovarianceMatrix(const CVData: TCMatrix; const WVar: TIVector);
var
  i, j, n: integer;
begin
  n := Length(CVData);
  writeln('Covariance Matrix   Variables: ', n);
  Write('                 ');
  for j := 0 to n - 1 do Write('Var', WVar[j]: 3, '     ');
  writeln;
  for i := 0 to n - 1 do begin
    Write('Var', WVar[i]: 3, '      ');
    for j := 0 to n - 1 do
      Write(CVData[i, j]: 11: 4);
    writeln;
  end;
  Pause;
  writeln;
end;

procedure DisplaySelectedData(const Data: TMatrix; const Condensed: Boolean);
var
  VData: TMatrix;
  WVar: TIVector;
  nVar: Integer;
begin
  SelectColumnsA(Data, VData, WVar, nVar);
  if not Condensed then
    DisplayData(VData, WVar)
  else
    DisplayDataCondensed(VData, WVar);
end;

procedure DisplayData(const PData: TMatrix; const PLabel: TIVector);
var
  i, j, m, n: Integer;
begin
  n := Length(PData);
  m := Length(PData[0].Value);
  writeln('Data     Obs: ', m, '   Vars: ', n);
  Write('      ');
  for i := 0 to n - 1 do write('     Col ', PLabel[i]: 3);
  writeln;
  for j := 0 to m - 1 do begin
    Write('Row', j + 1: 3);
    for i := 0 to n - 1 do
      Write(PData[i].Value[j]: 12: 5);
    writeln;
  end;
  Pause;
  writeln;
end;

procedure DisplayDataA(const PData: TMatrix; WVar: TString);
var
  i, j, m, n: Integer;
begin
  n := Length(PData);
  m := Length(PData[0].Value);
  writeln('Data     Obs: ', m, '   Vars: ', n);
  Write('  Var');
  for i := 0 to n - 1 do write('    ', PData[i].Name : 8);
  writeln;
  for j := 0 to m - 1 do begin
    Write('Obs', j + 1: 3);
    for i := 0 to n - 1 do
      Write(PData[i].Value[j]: 12: 5);
    writeln;
  end;
  Pause;
  writeln;
end;

procedure DisplayDataCondensed(const PData: TMatrix; const PLabel: TIVector);
var
  i, j, m, n: Integer;
begin
  n := Length(PData);
  m := Length(PData[0].Value);
  writeln('Data     Obs: ', m, '   Vars: ', n);
  Write('       ');
  for i := 0 to n - 1 do Write('C', PLabel[i]: 2, '    ');
  writeln;
  for j := 0 to m - 1 do begin
    Write('R', j + 1: 3);
    for i := 0 to n - 1 do
      Write('  ', PData[i].Value[j]: 5: 2);
    writeln;
  end;
  Pause;
  writeln;
end;

procedure DisplayMatrixM(const Mess: String; const PData: TMatrix);
var
  i, j, m, n: Integer;
begin
  n := Length(PData);
  m := Length(PData[0].Value);
  writeln(Mess, '     Rows: ', m, ' Columns: ', n);
  write('       ');
  for i := 0 to n - 1 do Write('C', i + 1: 2, '    ');
  writeln;
  for j := 0 to m - 1 do begin
    write('R', j + 1: 3);
    for i := 0 to n - 1 do
      write('  ', PData[i].Value[j]: 5: 2);
    writeln;
  end;
  Pause;
  writeln;
end;

procedure DME(const Mess: String; const Data: TMatrix);
var
  i, j, m, n: Integer;
begin
  n := Length(Data);
  m := Length(Data[0].Value);
  writeln(Mess, '     Rows: ', m, ' Columns: ', n);
  write('    ');
  for i := 0 to n - 1 do Write('C', i + 1: 2, '         ');
  writeln;
  for j := 0 to m - 1 do begin
    write('R', j + 1: 3);
    for i := 0 to n - 1 do
      write('  ', Data[i].Value[j]: 9: 5);
    writeln;
  end;
  Pause;
  writeln;
end;

procedure DE(const Data: TMatrix);
var
  i, j, m, n: Integer;
begin
  n := Length(Data);
  m := Length(Data[0].Value);
  writeln('         Rows: ', m, ' Columns: ', n);
  write('    ');
  for i := 0 to n - 1 do Write('C', i + 1: 2, '         ');
  writeln;
  for j := 0 to m - 1 do begin
    write('R', j + 1: 3);
    for i := 0 to n - 1 do
      write('  ', Data[i].Value[j]: 9: 5);
    writeln;
  end;
  Pause;
  writeln;
end;

procedure DisplayMatrix(const PData: TMatrix);
var                                             //add wvar
  i, j, m, n: Integer;
begin
  n := Length(PData);
  m := Length(PData[0].Value);
  writeln('Matrix     Rows: ', m, ' Columns: ', n);
  write('       ');
  for i := 0 to n - 1 do Write('C', i + 1: 2, '    ');
  writeln;
  for j := 0 to m - 1 do begin
    write('R', j + 1: 3);
    for i := 0 to n - 1 do
      write('  ', PData[i].Value[j]: 5: 2);
    writeln;
  end;
  Pause;
  writeln;
end;

procedure DisplayMatrixAll;
var
  i, j, m, n: Integer;
begin
  writeln('All Data  Rows: ', nRow, ' Columns: ', nCol);
  write('       ');
  for i := 0 to nCol - 1 do Write('C', i + 1: 2, '    ');
  writeln;
  for j := 0 to nRow - 1 do begin
    write('R', j + 1: 3);
    for i := 0 to nCol - 1 do
      write('  ', WData[i].Value[j]: 5: 2);
    writeln;
  end;
  Pause;
  writeln;
end;

procedure DC(const WVar: TIVector);
var
  j, m: Integer;
  Col: TVector;
begin
  m := Length(WVar);
  for j := 0 to m -1 do Col.Value[j] := WData[j].Value[WVar[0]];
  writeln(WData[WVar[0]].Name);
  for j := 0 to m - 1 do
    write('  ', Col.Value[j]: 5: 2);
  writeln;
  Pause;
end;

procedure DisplayColumn(const PData: TVector);
var                                             //add wvar
  j, n: Integer;
begin
  n := Length(PData.Value);
  writeln('Column     Rows: ', n);
  for j := 0 to n - 1 do
    write('  ', PData.Value[j]: 5: 2);
  writeln;
  Pause;
end;

procedure DisplayColumnM(const Mess: String; const PData: TVector);
var                                             //add wvar
  j, n: Integer;
begin
  n := Length(PData.Value);
  writeln(Mess, '     Rows: ', n);
  for j := 0 to n - 1 do
    write('  ', PData.Value[j]: 5: 2);
  writeln;
  Pause;
end;

procedure SaveData(const Mess: String; const Data: TMatrix);
var
  FileName: String;
  f: TextFile;
  Qinput: Char;
  i, j: Integer;
begin
  write('Save ', Mess, ' Data as CSV file? y/n ');
  readln(QInput);
  if UpCase(QInput) = 'Y' then begin
    write('Enter filename to save (e.g. matrix.csv): ');
    readln(FileName);
    AssignFile(f, FileName);
    rewrite(f);
    for j := 0 to High(Data) do begin
      for i := 0 to High(Data[j].Value) do begin
        write(f, Data[j].Value[i]: 0: 6);
        if i < High(Data[i].Value) then write(f, ',');
      end;
      writeln(f);
    end;
    closefile(f);
    writeln(' Data saved as ', FileName);
  end
  else
    writeln('Data not saved.');
  Pause;
end;

procedure SelectColumnsA(const Data: TMatrix; var VData: TMatrix;
  var WVar: TIVector; var nVar: Integer);
var
  i, j, k, m, n: Integer;
  VInput: Integer;
  Str: String;
begin
  n := Length(Data);
  m := Length(Data[0].Value);
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
  writeln;
  write('Selected ', nVar, ' variables: ');
  for i := 0 to nVar - 1 do
    write(' ', WVar[i]);
  writeln;
  // Construct VData from selected columns in WVar
  SetLength(VData, nCol);
  for j := 0 to m - 1 do
    for i := 0 to nVar - 1 do
      VData[i].Value[j] := Data[j].Value[WVar[i] - 1];
end;

{procedure CreateVName(const WData.Name: TMatrix; var VName: TString);
var
  i, j, k, m, n: Integer;
  Str: String;
begin
  m := Length(Data);
  n := Length(Data[0]);
  nVar := 0;
  writeln('Enter variable names. Hit Enter to end. Enter ''a'' for all.');
  repeat

    //Get input
    Write('>');
    Readln(Str); //I need the ln on read
//    Val(Str, VInput, k);

    //Exit on <CR>; all vars have been entered
    if Str = '' then
      if nVar = 0 then
        writeln('Input required.')
      else Break;

    //If 'A' then load all vars and exit
    if UpCase(Str) = 'A' then begin
      nVar := n;  //nVar is number of variables, up to n - 1
      SetLength(WData.Name, nVar);
      for k := 0 to nVar - 1 do //the Vars are 1..nVar
        VName[k] := WData.Name[k];
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
  writeln;
  write('Selected ', nVar, ' variables: ');
  for i := 0 to nVar - 1 do
    write(' ', WVar[i]);
  writeln;
  // Construct VData from selected columns in WVar
  SetLength(VData, m, nVar);
  for j := 0 to m - 1 do
    for i := 0 to nVar - 1 do
      VData[i, j] := Data[j, WVar[i] - 1];
end;}

procedure Trim(var Data: TMatrix; const iCol: Integer; const L, R: TFloat);
var                          // seems to work right now
  i, j, k: Integer;
begin
  {j := 0;
  repeat
    writeln('Data trim'); displaymatrix(Data);
    writeln('trim icol', icol - 1, ' j:', j); readln;
    writeln('on j ', j, '  r is ', r, 'Data is', Data[j, icol - 1]); }
//    if Data[iCol - 1].Value[j] > R then begin
      //writeln('in the k loop');
{      if j < High(Data) then begin
        for k := j to High(Data) - 1 do
          for i := 0 to High(Data[k]) do
            Data[k, i] := Data[k + 1, i];
      end;
      SetLength(Data, Length(Data) - 1);
    end
    else if Data[j, iCol - 1] < L then begin
      //writeln('in the k loop');
      if j < High(Data) then begin
        for k := j to High(Data) - 1 do
          for i := 0 to High(Data[k]) do
            Data[k, i] := Data[k + 1, i];
      end;
      SetLength(Data, Length(Data) - 1);
    end
    else Inc(j);
  until j > High(Data);
  if Length(Data) < 1 then
    writeln('Error on trimming: too few rows left.');}
end;

procedure Winsorize(var Data: TMatrix; iCol: Integer; lowerPct, upperPct: TFloat);
var                     // winsorizes column iCol - 1
  LowerBound, UpperBound: TFloat;
  j, n, lowerIdx, upperIdx: Integer;
  Col: TVector;
begin
  n := Length(Data);
  Col := Data[iCol - 1];
  //Sort
  QuickSort(Col, 0, High(Col.Value));
  // Step 2: Determine bounds
  lowerIdx := Round(lowerPct * n);
  upperIdx := Round(upperPct * n) - 1;
  lowerIdx := Max(0, Min(n - 1, lowerIdx));
  upperIdx := Max(0, Min(n - 1, upperIdx));
  lowerBound := Col.Value[lowerIdx];
  upperBound := Col.Value[upperIdx];
  // Step 3: Winsorize original Data
  for j := 0 to n - 1 do begin
    if Col.Value[j] < lowerBound then
      Col.Value[j] := lowerBound
    else if Col.Value[j] > upperBound then
      Col.Value[j] := upperBound;
  end;
  CheckWrite(Col); readln;
  //PutColumn(Data, iCol - 1, Col);
end;

procedure YeoJohnsonTransform(var Data: TMatrix; const iCol: Integer; Lambda: TFloat);
var
  j, m: Integer;
  xi: TFloat;
begin
  m := Length(Data);
  for j := 0 to m - 1 do begin
    xi := Data[iCol - 1].Value[j];
    if xi >= 0 then begin
      if Lambda = 0 then
        Data[iCol - 1].Value[j] := Ln(xi + 1)
      else
        Data[iCol - 1].Value[j] := (Power(xi + 1, Lambda) - 1) / Lambda;
    end
    else begin
      if Lambda = 2 then
        Data[iCol - 1].Value[j] := -Ln(-xi + 1)
      else
        Data[iCol - 1].Value[j] := - (Power(-xi + 1, 2 - Lambda) - 1) / (2 - Lambda);
    end;
  end;
end;

procedure ManuallyInputData(var WData: TMatrix; var nRow, nCol: Integer);
var
  i, j: Integer;
begin
  write('Number of columns (variables): ');
  readln(nCol);
  write('Number of rows (observations): ');
  readln(nRow);
  SetLength(WData, nCol);
  writeln('Enter each piece of Data, follow by <CR>.');
  for j := 0 to nRow - 1 do begin
    SetLength(WData[i].Value, nRow);
    write('Row ', j + 1);
    for i := 0 to nCol - 1 do begin
      write(' Col ', i + 1, ': ');
      read(WData[i].Value[j]);
//      write(' ');
    end;
//    writeln;
  end;
  writeln(' There are ', nRow, ' rows and ', nCol , nCol, ' columns.');
  writeln('Data were successfully entered.');
  Pause;
end;

procedure GetTestData(const Str: String);
var
  i: Integer;
begin
  Case Str of
    'T1': begin
      nRow := 12;
      nCol := 5;
      {SetLength(WData, nRow, nCol);
      SetLength(WData.Name, nCol);
      WData[0,0] := 2;   WData[0,1] := 3;   WData[0,2] := 4;   WData[0,3] := 5;   WData[0,4] := 6;
      WData[1,0] := 5;   WData[1,1] := 6;   WData[1,2] := 9;   WData[1,3] := 1;   WData[1,4] := 4;
      WData[2,0] := 6;   WData[2,1] := 4;   WData[2,2] := 3;   WData[2,3] := 11;  WData[2,4] := 14;
      WData[3,0] := 1;   WData[3,1] := 6;   WData[3,2] := 5;   WData[3,3] := 0;   WData[3,4] := 3;
      WData[4,0] := 1;   WData[4,1] := 5;   WData[4,2] := 3;   WData[4,3] := 14;  WData[4,4] := 17;
      WData[5,0] := 9;   WData[5,1] := 3;   WData[5,2] := 8;   WData[5,3] := 9;   WData[5,4] := 12;
      WData[6,0] := 9;   WData[6,1] := 3;   WData[6,2] := 16;  WData[6,3] := 0;   WData[6,4] := 3;
      WData[7,0] := 1;   WData[7,1] := 4;   WData[7,2] := 9;   WData[7,3] := 1;   WData[7,4] := 5;
      WData[8,0] := 10;  WData[8,1] := 5;   WData[8,2] := -3;  WData[8,3] := 13;  WData[8,4] := 16;
      WData[9,0] := 9;   WData[9,1] := 13;  WData[9,2] := 18;  WData[9,3] := 9;   WData[9,4] := 11;
      WData[10,0] := 9;  WData[10,1] := 3;  WData[10,2] := 16; WData[10,3] := 0;  WData[10,4] := 4;
      WData[11,0] := 11; WData[11,1] := 4;  WData[11,2] := -9; WData[11,3] := 1;  WData[11,4] := 5;
      WData.Name[0] := 'Ht'; WData.Name[1] := 'Wt'; WData.Name[2] := 'Size'; WData.Name[3] := 'A'; WData.Name[4] := 'b';
      writeln('Data created successfully.');
      Pause;                                }
    end;
    'T2': begin
      nRow := 2;
      nCol := 2;
      {SetLength(WData, nRow, nCol);   //invertible
      SetLength(WData.Name, nCol);
      WData.Name[0] := 'First'; WData.Name[1] := 'Second';
      WData[0,0] := 1; WData[0,1] := 0;
      WData[1,0] := 0; WData[1,1] := 2;
      writeln('Data created successfully.'); }
    end;
    'T3': begin
      nRow := 3;
      nCol := 3;
      SetLength(WData, nCol);             //invertible det =-9
      SetLength(WData[0].Name, 2);
      for i := 0 to nCol do SetLength(WData[i].Value, nRow);
      WData[0].Name := 'Way'; WData[1].Name := 'Here'; WData[2].Name := 'Avenue';
      WData[0].Value[0] := 2; WData[0].Value[1] := 1; WData[0].Value[2] := 3;
      WData[1].Value[0] := 0; WData[1].Value[1] := 1; WData[1].Value[2] := 4;
      WData[2].Value[0] := 5; WData[2].Value[1] := 2; WData[2].Value[2] := 1;
      writeln('Data created successfully.');
    end;
    'T4': begin
{      nRow := 5;
      nCol := 3;
      SetLength(WData, nRow, nCol);
      SetLength(WData.Name, nCol);
      WData.Name[0] := 'St'; WData.Name[1] := 'Yeah'; WData.Name[2] := 'Avenue';
      WData[0,0] := 2; WData[0,1] := 1; WData[0,2] := 3;
      WData[1,0] := 0; WData[1,1] := 1; WData[1,2] := 4;
      WData[2,0] := 9; WData[2,1] := 2; WData[2,2] := 1;
      WData[3,0] := 2; WData[3,1] := 7; WData[3,2] := 2;
      WData[4,0] := 1; WData[4,1] := 6; WData[4,2] := 1; }
      writeln('Data created successfully.');
    end;
  end;
end;

procedure ModifyData(var WData: TMatrix);
var
  WInput: Char;
  i, j, m, n, iCol, MultFactor, AddFactor: Integer;
  TrimR, TrimL, Lambda: TFloat;
  NormedData: TMatrix;
begin
  m := Length(WData);
  m := Length(WData[0].Value);
  repeat
    writeln('Options: ');
    writeln('1  Trim Data Based on One Variable (Remove Observations)');
    writeln('2  Normallize All Data (Mean=0 and SD=1)');
    writeln('3  Winsorize One Variable (Truncate Observations)');
    writeln('4  Yeo-Johnson Transform of One Variable (Modify Observations)');
    writeln('5  Create Linear Transform as New Variable)');
    writeln('X  Exit');
    write('>');
    readln(WInput);
    WInput := UpCase(WInput);
    Case WInput of
      '1': begin                       //no new Data matrix
        repeat
          write('Input number of variable: ');
          readln(iCol);      // iCol is the variable of interest
          if (iCol < 1) or (iCol > n) then writeln('Variable not in range.');
        until not((iCol < 1) or (iCol > n));
        write('Upper value: ');
        readln(TrimR);
        write('Lower value: ');
        readln(TrimL);
        Trim(WData, iCol, TrimL, TrimR);
        writeln('Data has been trimmed.');
      end;
{      '2': begin                        //new Data matrix created
        SetLength(NormedData, Length(WData), Length(WData[0]));
        NormalizeMatrix(WData, NormedData);
        for j := 0 to High(WData) do for i := 0 to High(WData[0]) do WData[i, j] := NormedData[i, j];
        writeln('Data has been normalized.');
      end;
      '3': begin                       //no new Data matrix
        repeat
          write('Input number of variable: ');
          readln(iCol);      // iCol is the variable of interest
          if (iCol < 1) or (iCol > n) then writeln('Variable not in range.');
        until not((iCol < 1) or (iCol > n));
        repeat
          write('Upper value between 0 and 1: ');
          readln(TrimR);
        until (TrimR > 0) and (TrimR < 0);
        repeat
          write('Lower value between 0 and 1: ');
          readln(TrimL);
        until (TrimR > 0) and (TrimR < 0);
        Winsorize(WData, iCol, TrimL, TrimR); //need to get Col, trim one column
        write('Variable ', iCol, ' has been winsorized.');
      end;
      '4': begin                       //no new Data matrix
        repeat
          write('Input number of variable: ');
          readln(iCol);      // iCol is the variable of interest
          if (iCol < 1) or (iCol > n) then writeln('Variable not in range.');
        until not((iCol < 1) or (iCol > n));
        write('Lambda: '); readln(Lambda);
        YeoJohnsonTransform(WData, iCol, Lambda);
        writeln('Data has been transformed.');
      end;
      '5': begin
        SetLength(WData, Length(WData), Length(WData[0]) + 1);
        write('Variable to weight: ');    //check that it's a variable
        readln(iCol);
        write('Factor to multiply Data by: ');
        readln(MultFactor);
        write('Factor to add to Data: ');
        readln(AddFactor);
        for j := 0 to m - 1 do
          WData[j, n + 1] := WData[j, iCol] * MultFactor + AddFactor;
        Inc(n);
      end;}
    end;
  until (WInput='X');
end;

end.
