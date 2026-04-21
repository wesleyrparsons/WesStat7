unit Parser;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Globals,
  Math,
  SysUtils;

function ParserExpression(const Pos: Integer; Line: String): WVectorType;

implementation

type
  TTokenType = (ttNone, ttIdent, ttNumber, ttPlus, ttMinus, ttMul, ttDiv, ttMod,
                ttPow, ttEq, ttNE, ttLT, ttLE, ttGT, ttGE,
                ttAnd, ttOr, ttNot, ttAbs, ttSqrt, ttRound,
                ttSort,
                ttCos, ttSin, ttATan, ttLog, ttLn,
                ttLParen, ttRParen, ttEnd);
  TToken = record
    Typ: TTokenType;
    Str: String;
    Val: FloatType;
  end;

  // See if a vector exists and return it.
function GetWVector(const Name: String; out Found: Boolean): WVectorType;
var
  i, j: Integer;
begin
  Found := False;
  SetLength(Result.Value, nRow);
  for i := 0 to nCol - 1 do                 // Check for variable.
    if UpCase(WData[i].Name) = UpCase(Name) then begin
      DeepCopyVector(WData[i], Result);     // Result := WData[i].
      Found := True;
      Break;
    end;

  if not Found then
    for i := 0 to nScal - 1 do              // Check for constant.
      if UpCase(SData[i].Name) = UpCase(Name) then begin
        for j := 0 to High(Result.Value) do
          Result.Value[j] := SData[i].Value;// Copy constant into GetWVector.
        Found := True;
        Break;
      end;
end;

// Create vector for intermediate resulsts and to return.
function CreateWVector(Val: FloatType): WVectorType; overload;
var
  i: Integer;
begin
  SetLength(Result.Value, nRow);
  for i := 0 to nRow - 1 do
    Result.Value[i] := Val;
  Result.Name := '';
end;

function CreateWVector(Arr: CVectorType): WVectorType; overload;
var i: Integer;
begin
  SetLength(Result.Value, Length(Arr));
  for i := 0 to High(Arr) do
    Result.Value[i] := Arr[i];
  Result.Name := '';
end;

// Parser function.
function ParserExpression(const Line: String): WVectorType;
var
  Token: TToken;

  // Get and consume the next token.
  procedure NextToken;
  var
    Ch: Char;
    S: String;
  begin
    Token.Typ := ttNone;
    Token.Str := '';
    Token.Val := 0.0;

    while (Pos <= Length(Line)) and (Line[Pos] in [' ', #9]) do Inc(Pos);
    if Pos > Length(Line) then begin Token.Typ := ttEnd; Exit; end;

    Ch := Line[Pos];
    if Ch in ['a'..'z','A'..'Z','_'] then begin
      S := '';
      // Build identifiers.
      while (Pos <= Length(Line)) and (Line[Pos] in ['a'..'z','A'..'Z','0'..'9','_']) do begin S := S + Line[Pos]; Inc(Pos); end;
      Token.Str := S;
      // Boolean and trig and other operators.
      if SameText(S, 'abs') then Token.Typ := ttAbs
      else if SameText(S, 'sqrt') then Token.Typ := ttSqrt
      else if SameText(S, 'round') then Token.Typ := ttRound
      else if SameText(S, 'cos') then Token.Typ := ttCos
      else if SameText(S, 'sin') then Token.Typ := ttSin
      else if SameText(S, 'atan') then Token.Typ := ttATan
      else if SameText(S, 'log') then Token.Typ := ttLog
      else if SameText(S, 'ln') then Token.Typ := ttLn
      else if SameText(S, 'and') then Token.Typ := ttAnd
      else if SameText(S, 'or') then Token.Typ := ttOr
      else if SameText(S, 'not') then Token.Typ := ttNot
      else if SameText(S, 'sort') then Token.Typ := ttSort
      else if SameText(S, 'mod') then Token.Typ := ttMod
      else Token.Typ := ttIdent;
    end
    // Numbers.
    else if Ch in ['0'..'9','.'] then begin
      S := '';
      while (Pos <= Length(Line)) and (Line[Pos] in ['0'..'9','.']) do
      begin S := S + Line[Pos]; Inc(Pos); end;
      Token.Typ := ttNumber;
      Token.Val := StrToFloat(S);
    end
    else begin
      Inc(Pos);
      // Arithmetic operators.
      case Ch of
        '+': Token.Typ := ttPlus;
        '-': Token.Typ := ttMinus;
        '*': Token.Typ := ttMul;
        '/': Token.Typ := ttDiv;
        '^': Token.Typ := ttPow;
        '(': Token.Typ := ttLParen;
        ')': Token.Typ := ttRParen;
        '=': Token.Typ := ttEq;
        '<': if (Pos <= Length(Line)) and (Line[Pos] = '=') then begin Inc(Pos); Token.Typ := ttLE; end
             else if (Pos <= Length(Line)) and (Line[Pos] = '>') then begin Inc(Pos); Token.Typ := ttNE; end
             else Token.Typ := ttLT;
        '>': if (Pos <= Length(Line)) and (Line[Pos] = '=') then begin Inc(Pos); Token.Typ := ttGE; end
             else Token.Typ := ttGT;
        else Token.Typ := ttNone;
      end;
    end;
  end;

  // Kinds of expressions.
  function Expression: WVectorType; forward;
  function Factor: WVectorType; forward;
  function Term: WVectorType; forward;
  function PowerExpr: WVectorType; forward;
  function Compare: WVectorType; forward;
  function Logic: WVectorType; forward;

  // Parse for negatives,
  function Unary: WVectorType;
  var
    V: WVectorType;
    i: Integer;
  begin
    if Token.Typ = ttMinus then begin
      NextToken;
      V := Factor;
      for i := 0 to High(V.Value) do
        V.Value[i] := -V.Value[i];
      Result := V;
    end
    else if Token.Typ = ttNot then begin
      NextToken;
      V := Factor;
      for i := 0 to High(V.Value) do
        V.Value[i] := Ord(V.Value[i] = 0.0);
      Result := V;
    end
    else
      Result := Factor;
  end;

  // Parse for numbers, identifiers, and parentheses, and operators.
  function Factor: WVectorType;
  var
    V: WVectorType;
    i: Integer;
    Found: Boolean;
    Func: TTokenType;
  begin
    if Token.Typ = ttNumber then begin                     // Numbers.
      Result := WVector_Create(Token.Val);
      NextToken;
    end
    else if Token.Typ = ttIdent then begin                 // Identifiers.
      Result := GetWVector(Token.Str, Found);
      if not Found then
        writeln('*', Token.Str, '* Error: variable not found in parse.');
      NextToken;
    end
    else if Token.Typ = ttLParen then begin                // Parentheses.
      NextToken;
      Result := Logic;
      if Token.Typ = ttRParen then NextToken;
    end
    else if Token.Typ = ttSort then begin                  // Sort.
      // New code added below for sort.
      NextToken;
      if Token.Typ <> ttLParen then Exit(CreateWVector(0));
      NextToken;
      V := Logic;
      if Token.Typ = ttRParen then NextToken;
      Quicksort(V, 0, High(V.Value));
      Result := V;
    end                                                    // Trig, etc.
    else if Token.Typ in [ttAbs, ttSqrt, ttRound, ttCos, ttSin, ttATan, ttLog, ttLn, ttSort] then begin
      Func := Token.Typ;
      NextToken;
      if Token.Typ <> ttLParen then Exit(CreateWVector(0));
      NextToken;
      V := Logic;
      if Token.Typ = ttRParen then NextToken;

      for i := 0 to High(V.Value) do
        case Func of
          ttAbs:   V.Value[i] := Abs(V.Value[i]);
          ttSqrt:  if V.Value[i] >= 0 then V.Value[i] := Sqrt(V.Value[i]) else V.Value[i] := NaN;
          ttRound: V.Value[i] := Round(V.Value[i]);
          ttCos:   V.Value[i] := Cos(V.Value[i]);
          ttSin:   V.Value[i] := Sin(V.Value[i]);
          ttATan:  V.Value[i] := ArcTan(V.Value[i]);
          ttLog:   if V.Value[i] > 0 then V.Value[i] := Log10(V.Value[i]) else V.Value[i] := NaN;
          ttLn:    if V.Value[i] > 0 then V.Value[i] := Ln(V.Value[i]) else V.Value[i] := NaN;
        end;
      Result := V;
    end
    else
      Result := CreateWVector(0);
  end;

  // Function for Multiplicatiom, Division, Mod.
  function Term: WVectorType;
  var
    W: WVectorType;
    i: Integer;
    Op: TTokenType;
  begin
    Result := Unary;
    while Token.Typ in [ttMul, ttDiv, ttMod] do begin
      Op := Token.Typ;
      NextToken;
      W := Unary;
      if Length(Result.Value) <> Length(W.Value) then Exit(CreateWVector(0));
      for i := 0 to High(Result.Value) do
        Case Op of
          // New Code below for mod.
          ttMul: Result.Value[i] := Result.Value[i] * W.Value[i];
          ttMod: Result.Value[i] := Result.Value[i] mod W.Value[i];
          ttDiv: if W.Value[i] <> 0 then
            Result.Value[i] := Result.Value[i] / W.Value[i];
        else
          Result.Value[i] := NaN;
        end;
    end;
  end;

  // Function for Power.
  function PowerExpr: WVectorType;
  var
    W: WVectorType;
    i: Integer;
  begin
    Result := Term;
    while Token.Typ = ttPow do begin
      NextToken;
      W := Term;
      if Length(Result.Value) <> Length(W.Value) then Exit(CreateWVector(0));

      for i := 0 to High(Result.Value) do
        Result.Value[i] := Power(Result.Value[i], W.Value[i]);
    end;
  end;

  // Function for Plus and Minus.
  function Expression: WVectorType;
  var
    W: WVectorType;
    i: Integer;
    Op: TTokenType;
  begin
    Result := PowerExpr;
    while Token.Typ in [ttPlus, ttMinus] do begin
      Op := Token.Typ;
      NextToken;
      W := PowerExpr;
      if Length(Result.Value) <> Length(W.Value) then Exit(CreateWVector(0));

      for i := 0 to High(Result.Value) do
        if Op = ttPlus then
          Result.Value[i] := Result.Value[i] + W.Value[i]
        else
          Result.Value[i] := Result.Value[i] - W.Value[i];
    end;
  end;

  // Function for Compares.
  function Compare: WVectorType;
  var
    W: WVectorType;
    i: Integer;
    Op: TTokenType;
  begin
    Result := Expression;
    if Token.Typ in [ttEq, ttNE, ttLT, ttLE, ttGT, ttGE] then begin
      Op := Token.Typ;
      NextToken;
      W := Expression;
      if Length(Result.Value) <> Length(W.Value) then Exit(CreateWVector(0));

      SetLength(Result.Value, Length(Result.Value));
      for i := 0 to High(Result.Value) do
        case Op of
          ttEq: Result.Value[i] := Ord(Abs(Result.Value[i] - W.Value[i]) < 1e-10);
          ttNE: Result.Value[i] := Ord(Abs(Result.Value[i] - W.Value[i]) >= 1e-10);
          ttLT: Result.Value[i] := Ord(Result.Value[i] < W.Value[i]);
          ttLE: Result.Value[i] := Ord(Result.Value[i] <= W.Value[i]);
          ttGT: Result.Value[i] := Ord(Result.Value[i] > W.Value[i]);
          ttGE: Result.Value[i] := Ord(Result.Value[i] >= W.Value[i]);
        end;
    end;
  end;

  // Function for Boolean.
  function Logic: WVectorType;
  var
    W: WVectorType;
    i: Integer;
    Op: TTokenType;
  begin
    Result := Compare;
    while Token.Typ in [ttAnd, ttOr] do begin
      Op := Token.Typ;
      NextToken;
      W := Compare;
      if Length(Result.Value) <> Length(W.Value) then Exit(CreateWVector(0));

      for i := 0 to High(Result.Value) do
        if Op = ttAnd then
          Result.Value[i] := Ord((Result.Value[i] <> 0) and (W.Value[i] <> 0))
        else
          Result.Value[i] := Ord((Result.Value[i] <> 0) or (W.Value[i] <> 0));
    end;
  end;

var
  VarName: string;
  Value: WVectorType;
begin
  if DebugOn then
    writeln('Entering ParserExpression. Line = *', Line, '*');
  Pos := 1;
  NextToken;

  if Token.Typ <> ttIdent then Exit(CreateWVector(0));
  VarName := Token.Str;
  NextToken;

  if Token.Typ <> ttEq then Exit(CreateWVector(0));
  NextToken;

  Value := Logic;

  if Token.Typ <> ttEnd then Exit(CreateWVector(0));
  DeepCopyVector(Value, Result);
  CleanUpWVector(Result);
  Result.Name := VarName;
end;

end.
