unit DataDisplay;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  Crt,
  Globals,
  Math,
  SysUtils;

{ Existence of Data }
procedure ReportDataState;
function NoDataExists: Boolean;

{ Number and Text Display }
function SmartFloat(const Value: Extended): String;
procedure WriteLeft(const s: String; const width: Integer);
function WriteTab(const Txt: String; const Value: FloatType): String;
function WriteTabln(const Txt: String; const Value: FloatType): String;
function CardinalToOrdinal(const n: Integer): String;
function CenterText(const S: string; Width: Integer): String;

{ Display Matrices }
procedure DisplayMatrixM(const Mess: String; const PData: WMatrixType; DisplayFormat: DisplayFormatType); overload;
procedure DisplayMatrixM(const Mess: String; const PData: CMatrixType; DisplayFormat: DisplayFormatType); overload;
procedure DisplayCorrelationMatrix(const Corr: CorrType; CData, SData: CMatrixType; const WVar: IVectorType);
procedure DisplayCovarianceMatrix(const CVData: CMatrixType; const WVar: IVectorType);

{ Display Vectors }
procedure DisplayWVar;
procedure DisplayVectorM(const Mess: String; const Col: CVectorType; DisplayFormat: DisplayFormatType); overload;
procedure DisplayVectorM(const Mess: String; const Col: WVectorType; DisplayFormat: DisplayFormatType); overload;

{ Display Scalars }
procedure DisplayScalarM(const Mess: String; const iScal: WScalarType;
  DisplayFormat: DisplayFormatType);
procedure DisplayAllScalars(DisplayFormat: DisplayFormatType);

implementation

var
  WVar: IVectorType;

procedure ReportDataState;
begin
  if not NoDataExists then
    writeln('There are  ', nRow, ' rows/observations, ', nCol, ' variables, and ', nScal, ' constants.');
end;

function NoDataExists: Boolean;
begin
  Result := False;
  if (WData = nil) and (SData = nil) then begin
    writeln('No data exist.');
    Result := True;
  end;
end;

{ Use SmartFloat in general number display }
{ Number and Text Display }
function SmartFloat(const Value: FloatType): String;
var
  v: FloatType;
begin
  // Handle special values
  if IsNan(Value) then Exit('NaN');
  if IsInfinite(Value) then begin
    if Value > 0 then
      Exit('Infinity')
    else
      Exit('-Infinity');
  end;

  // Normalize negative zero to plain 0
  v := Value;
  if (v = 0.0) then Exit('0');

  // Choose format based on magnitude
  if (Abs(v) < SmallThreshold) or (Abs(v) > LargeThreshold) then
    // Using Width and Precision, can be set by User
    // Use FloatToStrF with ffExponent for consistent exponent formatting
    Result := FloatToStrF(v, ffExponent, Width, Precision)
    // Example output: 1.23456E+03
  else begin
    // Fixed format
    Result := FloatToStrF(v, ffFixed, Width - 1, Precision);

    // Remove trailing zeros and possible trailing decimal point
    while (Length(Result) > 1) and (Result[Length(Result)] = '0') do
      Delete(Result, Length(Result), 1);
    if (Length(Result) > 0) and (Result[Length(Result)] = '.') then
      Delete(Result, Length(Result), 1);
  end;
end;

{ Use PSmartFloat in the display matrix procs }
function PSmartFloat(const Value: FloatType; const DisplayFormat: DisplayFormatType): String;
const
  SmallThreshold = 1e-6;
  LargeThreshold = 1e7;
begin

  // Handle special values
  if IsNan(Value) then Exit('NaN');
  if IsInfinite(Value) then begin
    if Value > 0 then
      Exit('Infinity')
    else
      Exit('-Infinity');
  end;

  // Normalize negative zero to plain 0
  if (Value = 0.0) then Exit('0');

  // Choose format based on magnitude
  Case DisplayFormat of
    T1: if (Abs(Value) < SmallThreshold) or (Abs(Value) > LargeThreshold) then
          Result := Format('%8.3e', [Value])
        else begin
          // Fixed format: up to 3 decimals, remove trailing zeros
          Result := FloatToStrF(Value, ffFixed, 8, 3);
          // Remove trailing zeros and possible trailing decimal point
          while (Length(Result) > 1) and (Result[Length(Result)] = '0') do
            Delete(Result, Length(Result), 1);
          if (Length(Result) > 0) and (Result[Length(Result)] = '.') then
            Delete(Result, Length(Result), 1);
        end;
    T2: if (Abs(Value) < SmallThreshold) or (Abs(Value) > LargeThreshold) then
          Result := Format('%11.4e', [Value])
        else begin
          // Fixed format: up to 4 decimals, remove trailing zeros
          Result := FloatToStrF(Value, ffFixed, 11, 4);
          // Remove trailing zeros and possible trailing decimal point
          while (Length(Result) > 1) and (Result[Length(Result)] = '0') do
            Delete(Result, Length(Result), 1);
          if (Length(Result) > 0) and (Result[Length(Result)] = '.') then
            Delete(Result, Length(Result), 1);
        end;
    T3: if (Abs(Value) < SmallThreshold) or (Abs(Value) > LargeThreshold) then
          Result := Format('%15.6e', [Value])
        else begin
          // Fixed format: up to 6 decimals, remove trailing zeros
          Result := FloatToStrF(Value, ffFixed, 15, 6);
          // Remove trailing zeros and possible trailing decimal point
          while (Length(Result) > 1) and (Result[Length(Result)] = '0') do
            Delete(Result, Length(Result), 1);
          if (Length(Result) > 0) and (Result[Length(Result)] = '.') then
            Delete(Result, Length(Result), 1);
        end;
  end;
end;

procedure WriteLeft(const s: String; const width: Integer);
begin
  Write(Copy(s + StringOfChar(' ', width), 1, width));
end;

function WriteTabln(const Txt: String; const Value: FloatType): String;
var
  Num, Padded: String;
  PadLen: Integer;
  RightAlign: Boolean = False;
  FieldWidth: Integer = 0;
begin
  Num := SmartFloat(Value);
  Result := Txt + Num;

  if FieldWidth > 0 then begin
    PadLen := FieldWidth - Length(Num);
    if PadLen <= 0 then
      Padded := Num
    else if RightAlign then
      Padded := StringOfChar(' ', PadLen) + Num
    else
      Padded := Num + StringOfChar(' ', PadLen);
    Result := Txt + Padded;
  end;
  writeln(Result);
end;

function WriteTab(const Txt: String; const Value: FloatType): String;
var
  Num, Padded: String;
  Pad, PadLen: Integer;
  RightAlign: Boolean = False;
  FieldWidth: Integer = 0;
begin
  Num := SmartFloat(Value);
  Result := Txt + Num;

  if FieldWidth > 0 then begin
    PadLen := FieldWidth - Length(Num);
    if PadLen <= 0 then
      Padded := Num
    else if RightAlign then
      Padded := StringOfChar(' ', PadLen) + Num
    else
      Padded := Num + StringOfChar(' ', PadLen);
    Result := Txt + Padded;
  end;

  Pad := Length(Num) + Length(Txt);
  if Pad > (Tab - 1) then
    Pad := Tab
  else
    Result := Txt + Num + StringofChar(' ', Tab - Pad);
  Write(Result);
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

function CenterText(const S: string; Width: Integer): String;
var
  L, PadTotal, PadLeft, PadRight: Integer;
begin
  if Width <= 0 then Exit('');

  L := Length(S);
  if L >= Width then
    // Trim to width (preserve leftmost characters).
    Result := Copy(S, 1, Width)
  else begin
    PadTotal := Width - L;
    PadLeft := PadTotal div 2;
    PadRight := PadTotal - PadLeft;
    Result := StringOfChar(' ', PadLeft) + S + StringOfChar(' ', PadRight);
  end;
end;

procedure DisplayMatrixM(const Mess: String; const PData: WMatrixType; DisplayFormat: DisplayFormatType); overload;
var
  i, j, m, n, v: Integer;
  s: String;
begin
  if NoDataExists then Exit;
  n := Length(PData);
  m := Length(PData[0].Value);
  v := ord(DisplayFormat)*3 + 9;

  writeln(Mess, '     Variables     Rows: ', m, '    Columns: ', n);
  write(' ' : 9);
  for i := 0 to n - 1 do begin
    // Ensure Name fits in v - 1 characters.
    s := PData[i].Name;
    if Length(s) > v - 2 then
      s := Copy(s, 1, v - 2);
    Write(' ', s: v);
  end;
  writeln;

  for j := 0 to m - 1 do begin
    write('Row ', (j + 1) : 4, ' ');
    for i := 0 to n - 1 do
      write(' ', PSmartFloat(PData[i].Value[j], DisplayFormat): v);
    writeln;
  end;
end;

procedure DisplayMatrixM(const Mess: String; const PData: CMatrixType; DisplayFormat: DisplayFormatType); overload;
var
  i, j, m, n, v: Integer;
  s: String;
begin
  if NoDataExists then Exit;
  n := Length(PData);
  m := Length(PData[0]);
  v := Ord(T1) * 3 + 9;

  writeln(Mess, '     Rows: ', m, '    Columns: ', n);
  write(' ' : 9);
  for i := 0 to n - 1 do begin
    // Ensure Name fits in v - 1 characters.
    s := 'Col' + IntToStr(i + 1);
    if Length(s) > v - 2 then
      s := Copy(s, 1, v - 2);
    Write(' ', s: v);
  end;
  writeln;

  for j := 0 to m - 1 do begin
    write('Row ', (j + 1) : 4, ' ');
    for i := 0 to n - 1 do
      write(' ', PSmartFloat(PData[i, j], DisplayFormat): v);
    writeln;
  end;
end;

procedure DisplayCorrelationMatrix(const Corr: CorrType; CData, SData: CMatrixType; const WVar: IVectorType);
var
  i, j, n: integer;
begin
  if NoDataExists then Exit;
  n := Length(CData);
  Case Corr of
    Pearson: write('Pearson''s R');
    Spearman: write('Spearman''s Rho');
    Kendall: write('Kendall''s Tau');
    Hoeffding: write('Hoeffding''s D');
    MI: write('Mutual Information');
    NMI: write('Normalized Mutual Information');
  end;

  writeln(' Correlation Matrix       Variables: ', n);
  write('          ');
  for j := 0 to n - 1 do
    Write(WData[WVar[j]].Name : 10);
  writeln;

  for i := 0 to n - 1 do begin
    WriteLeft(WData[WVar[i]].Name, 10);
    for j := 0 to n - 1 do
      Write(CData[i, j]: 10: Precision);
    writeln;
    if ShowSignificance and not (Corr = MI) and not (Corr = NMI) then begin
      Write('          ');
      for j := 0 to n - 1 do
        Write('  (', SData[i, j] : 6: Precision - 1, ')');
      writeln;
    end;
  end;

  if ShowSignificance and not (Corr = MI) and not (Corr = NMI) then
    writeln('Two-tailed significance in parentheses.');
end;

procedure DisplayCovarianceMatrix(const CVData: CMatrixType; const WVar: IVectorType);
var
  i, j, n: integer;
begin
  if NoDataExists then Exit;
  n := Length(CVData);
  writeln('Covariance Matrix       Variables: ', n);
  write('          ');
  for j := 0 to n - 1 do
    Write(WData[WVar[j]].Name : Width + 4);
  writeln;

  for i := 0 to n - 1 do begin
    WriteLeft(WData[WVar[i]].Name, 10);
    for j := 0 to n - 1 do
      Write(CVData[i, j] : Width + 4: Precision);
    writeln;
  end;
end;

{ Display Vectors }
procedure DisplayWVar;
var
  i: Integer;
begin
  write('Length = ', Length(WVar), '     WVar: ');
  for i := 0 to High(WVar)   do
    write(WVar[i]: 3);
  writeln;
end;

procedure DisplayVectorM(const Mess: String; const Col: WVectorType;
  DisplayFormat: DisplayFormatType); overload;
var
  j, m: Integer;
begin
  if NoDataExists then Exit;
  m := Length(Col.Value);
  writeln(Mess, '   Rows: ', m);
  if DisplayFormat in [T2, T3] then begin
    write('Mean = ', SmartFloat(Col.Mean), '  ');
    write('Median = ', SmartFloat(Col.Median), '  ');
    write('Std. Dev. = ', SmartFloat(Col.StdDev), '  ');
    write('Minimum = ', SmartFloat(Col.Min), '  ');
    write('Maximum = ', SmartFloat(Col.Max), '  ');
  end;
  write('Values: ');
  for j := 0 to m - 1 do
    Case DisplayFormat of
      T1 : write('  ', Col.Value[j]: 5: 2);
      T2 : write('  ', Col.Value[j]: 8: 5);
      T3 : write('  ', Col.Value[j]: Width: Precision);
    end;
  writeln;
end;

procedure DisplayVectorM(const Mess: String; const Col: CVectorType;
  DisplayFormat: DisplayFormatType); overload;
var
  j, m: Integer;
begin
  if NoDataExists then Exit;
  m := Length(Col);
  writeln(Mess, ' Rows: ', m);
  for j := 0 to m - 1 do
    Case DisplayFormat of
      T1 : write('  ', Col[j]: 5: 2);
      T2 : write('  ', Col[j]: 8: 5);
      T3 : write('  ', Col[j]: Width: Precision);
    end;
  writeln;
end;

{ Display Scalar }
procedure DisplayScalarM(const Mess: String; const iScal: WScalarType;
  DisplayFormat: DisplayFormatType);
begin
  if NoDataExists then Exit;
  write(Mess, iScal.Name, ': ', PSmartFloat(iScal.Value, DisplayFormat));
  writeln;
end;

procedure DisplayAllScalars(DisplayFormat: DisplayFormatType);
var
  i: Integer;
begin
  if NoDataExists then Exit;
  writeln('Constants: ', nScal);
  for i := 0 to nScal - 1 do
    write(SData[i].Name, ': ', PSmartFloat(SData[i].Value, DisplayFormat), '   ');
end;

end.
