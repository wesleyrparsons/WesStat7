unit EE;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, Math, SysUtils;

implementation

type
  TValueKind = (vkScalar, vkVector, vkMatrix, vkBoolean);

  TDoubleVector = array of Double;
  TDoubleMatrix = array of TDoubleVector;

  TValue = record
    kind: TValueKind;
    scalar: Double;
    bool: Boolean;
    vec: TDoubleVector;
    mat: TDoubleMatrix;
  end;

  TVarEntry = class
    Value: TValue;
  end;

  type
  TParserr = class
  public
    FExpr: string;
    FPos: Integer;
    function Peek: Char;
    function Next: Char;
    procedure SkipSpaces;
    procedure Expect(ch: Char);
    function ParserStatement: TValue;
    function ParserExpression: TValue;
    function ParserLogic: TValue;
    function ParserComparison: TValue;
    function ParserAddSub: TValue;
    function ParserTerm: TValue;
    function ParserFactor: TValue;
    function ParserUnary: TValue;
    function ParserPrimary: TValue;
    function ParserIdentifier: TValue;
    function ParserNumber: TValue;
    function ParserArrayLiteral: TValue; // vector or matrix
    function ParserFunction(const Name: string): TValue;
  public
    function ParserTop: TValue;
  end;

  { ---------- Parserr Implementation ---------- }

function TParserr.Peek: Char;
begin
  if FPos > Length(FExpr) then
    Exit(#0)
  else
    Exit(FExpr[FPos]);
end;

function TParserr.Next: Char;
begin
  if FPos < Length(FExpr) then
    Inc(FPos)
  else
    FPos := Length(FExpr)+1;
  Result := Peek;
end;

procedure TParserr.SkipSpaces;
begin
  while (Peek in [' ', #9, #10, #13]) do Next;
end;

procedure TParserr.Expect(ch: Char);
begin
  SkipSpaces;
  if Peek <> ch then raise Exception.CreateFmt('Expected "%s" at pos %d', [ch, FPos]);
  Next; // consume
end;

function MakeScalar(x: Double): TValue; forward;
function VarGet(const name: string; out found: Boolean): TValue;   forward;
function TParserr.ParserFunction(const Name: string): TValue; forward;

function TParserr.ParserNumber: TValue;
var start: Integer; s: string;
begin
  SkipSpaces;
  start := FPos;
  // simple numeric scanner: digits, dot, e/E and signs inside exponent
  while (Peek in ['0'..'9', '.', 'e', 'E', '+', '-']) do
  begin
    Inc(FPos);
    if FPos > Length(FExpr) then Break;
  end;
  // we moved one too far with loop; correct
  if FPos > Length(FExpr) then s := Copy(FExpr, start, Length(FExpr)-start+1)
  else s := Copy(FExpr, start, FPos - start);
  // set FPos to current char (we will not consume extra here)
  Result := MakeScalar(StrToFloatDef(s, 0.0));
end;

function TParserr.ParserIdentifier: TValue;
var name: string; found: Boolean;
begin
  SkipSpaces;
  name := '';
  while Peek in ['A'..'Z','a'..'z','0'..'9','_'] do
  begin
    name := name + Peek;
    Next;
  end;
  SkipSpaces;
  if Peek = '(' then
    Exit(ParserFunction(name));
  Result := VarGet(name, found);
end;

function TParserr.ParserArrayLiteral: TValue;
{ Accepts either:
    [ expr, expr, ... ]    -> vector
    [ [r1c1, r1c2], [r2c1, r2c2] ] -> matrix
  Rows must be consistent.
}


procedure EvalExpression(const Expr: string);


function EvalToValue(const Expr: string): TValue;


var
  Vars: TStringList;

{ ---------- Helpers / Constructors ---------- }

function MakeScalar(x: Double): TValue;
begin
  Result.kind := vkScalar;
  Result.scalar := x;
  Result.bool := False;
  SetLength(Result.vec, 0);
  SetLength(Result.mat, 0);
end;

function MakeBool(b: Boolean): TValue;
begin
  Result.kind := vkBoolean;
  Result.bool := b;
  Result.scalar := Ord(b);
  SetLength(Result.vec, 0);
  SetLength(Result.mat, 0);
end;

function MakeVectorFromArray(const a: TDoubleVector): TValue;
begin
  Result.kind := vkVector;
  Result.vec := Copy(a, 0, Length(a));
  Result.scalar := 0.0;
  Result.bool := False;
  SetLength(Result.mat, 0);
end;

function MakeMatrixFromArray(const m: TDoubleMatrix): TValue;
var i: Integer;
begin
  Result.kind := vkMatrix;
  SetLength(Result.mat, Length(m));
  for i := 0 to High(m) do
    Result.mat[i] := Copy(m[i], 0, Length(m[i]));
  Result.scalar := 0.0;
  Result.bool := False;
  SetLength(Result.vec, 0);
end;

{ ---------- Symbol table ---------- }

procedure VarSet(const name: string; const value: TValue);
var
  idx: Integer; entry: TVarEntry;
begin
  idx := Vars.IndexOf(name);
  if idx = -1 then
  begin
    entry := TVarEntry.Create;
    entry.Value := value;
    Vars.AddObject(name, entry);
  end
  else
    TVarEntry(Vars.Objects[idx]).Value := value;
end;

function VarGet(const name: string; out found: Boolean): TValue;
var idx: Integer;
begin
  idx := Vars.IndexOf(name);
  found := idx <> -1;
  if found then
    Result := TVarEntry(Vars.Objects[idx]).Value
  else
    raise Exception.Create('Unknown identifier: ' + name);
end;

{ ---------- elementwise utilities ---------- }

procedure CheckVecSizes(const A, B: TDoubleVector);
begin
  if Length(A) <> Length(B) then
    raise Exception.Create('Vector length mismatch');
end;

procedure CheckMatSizes(const A, B: TDoubleMatrix);
var i: Integer;
begin
  if Length(A) <> Length(B) then
    raise Exception.Create('Matrix row-count mismatch');
  if (Length(A) > 0) and (Length(B) > 0) and (Length(A[0]) <> Length(B[0])) then
    raise Exception.Create('Matrix column-count mismatch');
  // Also check each row length
  i := 0;
  for i := 0 to High(A) do
    if Length(A[i]) <> Length(A[0]) then
      raise Exception.Create('Irregular matrix (A) - rows unequal');
  for i := 0 to High(B) do
    if Length(B[i]) <> Length(B[0]) then
      raise Exception.Create('Irregular matrix (B) - rows unequal');
end;

function ElementwiseVecOp(const A, B: TDoubleVector; op: Char): TDoubleVector;
var i: Integer;
begin
  CheckVecSizes(A, B);
  SetLength(Result, Length(A));
  for i := 0 to High(A) do
    case op of
      '+': Result[i] := A[i] + B[i];
      '-': Result[i] := A[i] - B[i];
      '*': Result[i] := A[i] * B[i];
      '/': if B[i]=0 then begin raise Exception.Create('Division by zero in vector element'); Result[i] := A[i] / B[i]; end;
      '^': Result[i] := Power(A[i], B[i]);
      '%': Result[i] := Fmod(A[i], B[i]);
    else raise Exception.Create('Unknown elementwise op for vectors');
    end;
end;

function ElementwiseMatOp(const A, B: TDoubleMatrix; op: Char): TDoubleMatrix;
var i,j: Integer;
begin
  CheckMatSizes(A,B);
  SetLength(Result, Length(A));
  for i := 0 to High(A) do
  begin
    SetLength(Result[i], Length(A[i]));
    for j := 0 to High(A[i]) do
      case op of
        '+': Result[i][j] := A[i][j] + B[i][j];
        '-': Result[i][j] := A[i][j] - B[i][j];
        '*': Result[i][j] := A[i][j] * B[i][j];
        '/': if B[i][j]=0 then begin raise Exception.Create('Division by zero in matrix element'); Result[i][j] := A[i][j] / B[i][j]; end;
        '^': Result[i][j] := Power(A[i][j], B[i][j]);
        '%': Result[i][j] := Fmod(A[i][j], B[i][j]);
      else raise Exception.Create('Unknown elementwise op for matrices');
      end;
  end;
end;

function ScalarVecOp(s: Double; const V: TDoubleVector; op: Char; leftScalar: Boolean): TDoubleVector;
var i: Integer;
begin
  SetLength(Result, Length(V));
  for i := 0 to High(V) do
    case op of
      '+': if leftScalar then Result[i] := s + V[i] else Result[i] := V[i] + s;
      '-': if leftScalar then Result[i] := s - V[i] else Result[i] := V[i] - s;
      '*': if leftScalar then Result[i] := s * V[i] else Result[i] := V[i] * s;
      '/': if leftScalar then Result[i] := s / V[i] else Result[i] := V[i] / s;
      '^': if leftScalar then Result[i] := Power(s, V[i]) else Result[i] := Power(V[i], s);
      '%': if leftScalar then Result[i] := Fmod(s, V[i]) else Result[i] := Fmod(V[i], s);
    else raise Exception.Create('Unknown scalar-vector op');
    end;
end;

function ScalarMatOp(s: Double; const M: TDoubleMatrix; op: Char; leftScalar: Boolean): TDoubleMatrix;
var i,j: Integer;
begin
  SetLength(Result, Length(M));
  for i := 0 to High(M) do
  begin
    SetLength(Result[i], Length(M[i]));
    for j := 0 to High(M[i]) do
      case op of
        '+': if leftScalar then Result[i][j] := s + M[i][j] else Result[i][j] := M[i][j] + s;
        '-': if leftScalar then Result[i][j] := s - M[i][j] else Result[i][j] := M[i][j] - s;
        '*': if leftScalar then Result[i][j] := s * M[i][j] else Result[i][j] := M[i][j] * s;
        '/': if leftScalar then Result[i][j] := s / M[i][j] else Result[i][j] := M[i][j] / s;
        '^': if leftScalar then Result[i][j] := Power(s, M[i][j]) else Result[i][j] := Power(M[i][j], s);
        '%': if leftScalar then Result[i][j] := Fmod(s, M[i][j]) else Result[i][j] := Fmod(M[i][j], s);
      else raise Exception.Create('Unknown scalar-matrix op');
      end;
  end;
end;

{ Matrix multiply (classic A * B) - not elementwise }
function MatMul(const A, B: TDoubleMatrix): TDoubleMatrix;
var i,j,k, nArows, nAcols, nBrows, nBcols: Integer;
    sum: Double;
begin
  nArows := Length(A);
  nAcols := 0; if nArows>0 then nAcols := Length(A[0]);
  nBrows := Length(B);
  nBcols := 0; if nBrows>0 then nBcols := Length(B[0]);
  if nAcols <> nBrows then raise Exception.Create('Matrix multiply dimension mismatch');
  SetLength(Result, nArows);
  for i := 0 to nArows-1 do
  begin
    SetLength(Result[i], nBcols);
    for j := 0 to nBcols-1 do
    begin
      sum := 0.0;
      for k := 0 to nAcols-1 do
        sum := sum + A[i][k] * B[k][j];
      Result[i][j] := sum;
    end;
  end;
end;

{ ---------- Arithmetic dispatching ---------- }

function AddVal(const A, B: TValue): TValue;
begin
  if (A.kind = vkScalar) and (B.kind = vkScalar) then Exit(MakeScalar(A.scalar + B.scalar));
  if (A.kind = vkVector) and (B.kind = vkVector) then Exit(MakeVectorFromArray(ElementwiseVecOp(A.vec, B.vec, '+')));
  if (A.kind = vkMatrix) and (B.kind = vkMatrix) then Exit(MakeMatrixFromArray(ElementwiseMatOp(A.mat, B.mat, '+')));
  if (A.kind = vkScalar) and (B.kind = vkVector) then Exit(MakeVectorFromArray(ScalarVecOp(A.scalar, B.vec, '+', True)));
  if (A.kind = vkVector) and (B.kind = vkScalar) then Exit(MakeVectorFromArray(ScalarVecOp(B.scalar, A.vec, '+', False)));
  if (A.kind = vkScalar) and (B.kind = vkMatrix) then Exit(MakeMatrixFromArray(ScalarMatOp(A.scalar, B.mat, '+', True)));
  if (A.kind = vkMatrix) and (B.kind = vkScalar) then Exit(MakeMatrixFromArray(ScalarMatOp(B.scalar, A.mat, '+', False)));
  raise Exception.Create('Unsupported addition types');
end;

function SubVal(const A, B: TValue): TValue;
begin
  if (A.kind = vkScalar) and (B.kind = vkScalar) then Exit(MakeScalar(A.scalar - B.scalar));
  if (A.kind = vkVector) and (B.kind = vkVector) then Exit(MakeVectorFromArray(ElementwiseVecOp(A.vec, B.vec, '-')));
  if (A.kind = vkMatrix) and (B.kind = vkMatrix) then Exit(MakeMatrixFromArray(ElementwiseMatOp(A.mat, B.mat, '-')));
  if (A.kind = vkScalar) and (B.kind = vkVector) then Exit(MakeVectorFromArray(ScalarVecOp(A.scalar, B.vec, '-', True)));
  if (A.kind = vkVector) and (B.kind = vkScalar) then Exit(MakeVectorFromArray(ScalarVecOp(B.scalar, A.vec, '-', False)));
  if (A.kind = vkScalar) and (B.kind = vkMatrix) then Exit(MakeMatrixFromArray(ScalarMatOp(A.scalar, B.mat, '-', True)));
  if (A.kind = vkMatrix) and (B.kind = vkScalar) then Exit(MakeMatrixFromArray(ScalarMatOp(B.scalar, A.mat, '-', False)));
  raise Exception.Create('Unsupported subtraction types');
end;

function MulVal(const A, B: TValue): TValue;
var
  i, j: Integer;
  sum: Float;
  resvec: array of Double;
begin
  if (A.kind = vkScalar) and (B.kind = vkScalar) then Exit(MakeScalar(A.scalar * B.scalar));
  if (A.kind = vkVector) and (B.kind = vkVector) then Exit(MakeVectorFromArray(ElementwiseVecOp(A.vec, B.vec, '*')));
  if (A.kind = vkMatrix) and (B.kind = vkMatrix) then Exit(MakeMatrixFromArray(ElementwiseMatOp(A.mat, B.mat, '*')));
  if (A.kind = vkScalar) and (B.kind = vkVector) then Exit(MakeVectorFromArray(ScalarVecOp(A.scalar, B.vec, '*', True)));
  if (A.kind = vkVector) and (B.kind = vkScalar) then Exit(MakeVectorFromArray(ScalarVecOp(B.scalar, A.vec, '*', False)));
  if (A.kind = vkScalar) and (B.kind = vkMatrix) then Exit(MakeMatrixFromArray(ScalarMatOp(A.scalar, B.mat, '*', True)));
  if (A.kind = vkMatrix) and (B.kind = vkScalar) then Exit(MakeMatrixFromArray(ScalarMatOp(B.scalar, A.mat, '*', False)));

  // optionally allow matrix*vector (linear algebra)
  if (A.kind = vkMatrix) and (B.kind = vkVector) then
begin
  if (Length(A.mat) = 0) then raise Exception.Create('Matrix has no rows');
  if Length(A.mat[0]) <> Length(B.vec) then
    raise Exception.Create('Matrix/Vector dimension mismatch');
  SetLength(resVec, Length(A.mat));
  for i := 0 to High(A.mat) do
  begin
    sum := 0.0;
    for j := 0 to High(A.mat[i]) do
      sum := sum + A.mat[i][j] * B.vec[j];
    resVec[i] := sum;
  end;
  Exit(MakeVectorFromArray(resVec));
end;


  // allow matmul via function 'matmul' - not here
  raise Exception.Create('Unsupported multiplication types');
end;

function DivVal(const A, B: TValue): TValue;
begin
  if (A.kind = vkScalar) and (B.kind = vkScalar) then Exit(MakeScalar(A.scalar / B.scalar));
  if (A.kind = vkVector) and (B.kind = vkVector) then Exit(MakeVectorFromArray(ElementwiseVecOp(A.vec, B.vec, '/')));
  if (A.kind = vkMatrix) and (B.kind = vkMatrix) then Exit(MakeMatrixFromArray(ElementwiseMatOp(A.mat, B.mat, '/')));
  if (A.kind = vkScalar) and (B.kind = vkVector) then Exit(MakeVectorFromArray(ScalarVecOp(A.scalar, B.vec, '/', True)));
  if (A.kind = vkVector) and (B.kind = vkScalar) then Exit(MakeVectorFromArray(ScalarVecOp(B.scalar, A.vec, '/', False)));
  if (A.kind = vkScalar) and (B.kind = vkMatrix) then Exit(MakeMatrixFromArray(ScalarMatOp(A.scalar, B.mat, '/', True)));
  if (A.kind = vkMatrix) and (B.kind = vkScalar) then Exit(MakeMatrixFromArray(ScalarMatOp(B.scalar, A.mat, '/', False)));
  raise Exception.Create('Unsupported division types');
end;

function PowVal(const A, B: TValue): TValue;
begin
  if (A.kind = vkScalar) and (B.kind = vkScalar) then Exit(MakeScalar(Power(A.scalar, B.scalar)));
  if (A.kind = vkVector) and (B.kind = vkVector) then Exit(MakeVectorFromArray(ElementwiseVecOp(A.vec, B.vec, '^')));
  if (A.kind = vkMatrix) and (B.kind = vkMatrix) then Exit(MakeMatrixFromArray(ElementwiseMatOp(A.mat, B.mat, '^')));
  if (A.kind = vkScalar) and (B.kind = vkVector) then Exit(MakeVectorFromArray(ScalarVecOp(A.scalar, B.vec, '^', True)));
  if (A.kind = vkVector) and (B.kind = vkScalar) then Exit(MakeVectorFromArray(ScalarVecOp(B.scalar, A.vec, '^', False)));
  if (A.kind = vkScalar) and (B.kind = vkMatrix) then Exit(MakeMatrixFromArray(ScalarMatOp(A.scalar, B.mat, '^', True)));
  if (A.kind = vkMatrix) and (B.kind = vkScalar) then Exit(MakeMatrixFromArray(ScalarMatOp(B.scalar, A.mat, '^', False)));
  raise Exception.Create('Unsupported power types');
end;

{ ---------- Functions ---------- }

function DotProduct(const A, B: TValue): TValue;
var i: Integer; sum: Double;
begin
  if (A.kind <> vkVector) or (B.kind <> vkVector) then raise Exception.Create('dot requires two vectors');
  CheckVecSizes(A.vec, B.vec);
  sum := 0.0;
  for i := 0 to High(A.vec) do sum := sum + A.vec[i]*B.vec[i];
  Result := MakeScalar(sum);
end;

function MatMulFunction(const A, B: TValue): TValue;
var c: Integer;
begin
  if (A.kind <> vkMatrix) or (B.kind <> vkMatrix) then raise Exception.Create('matmul requires two matrices');
  Result := MakeMatrixFromArray(MatMul(A.mat, B.mat));
end;

var
  rows: TDoubleMatrix;
  rowVec: TDoubleVector;
  isMatrix: Boolean;
  firstElemKind: TValueKind;
  tmpVal: TValue;
  c, r: Integer;
begin
  SkipSpaces;
  Expect('[');
  SetLength(rows, 0);
  SetLength(rowVec, 0);
  isMatrix := False;
  firstElemKind := vkScalar;

  SkipSpaces;
  if Peek = ']' then
  begin
    // empty vector
    Expect(']');
    Exit(MakeVectorFromArray(rowVec));
  end;

  // Detect whether first element is '[' => matrix
  if Peek = '[' then
  begin
    isMatrix := True;
    // parse rows
    while True do
    begin
      // parse one row: [ expr, expr, ... ]
      Expect('[');
      SetLength(rowVec, 0);
      while True do
      begin
        tmpVal := ParserExpression;
        if tmpVal.kind <> vkScalar then
          raise Exception.Create('Matrix elements must be scalars');
        SetLength(rowVec, Length(rowVec)+1);
        rowVec[High(rowVec)] := tmpVal.scalar;
        SkipSpaces;
        if Peek = ',' then Next else Break;
      end;
      Expect(']');
      // append row
      SetLength(rows, Length(rows)+1);
      rows[High(rows)] := Copy(rowVec, 0, Length(rowVec));
      SkipSpaces;
      if Peek = ',' then Next else Break;
    end;
    Expect(']');
    // validate equal row lengths
    if Length(rows) > 0 then
    begin
      c := Length(rows[0]);
      for r := 0 to High(rows) do
        if Length(rows[r]) <> c then raise Exception.Create('Matrix rows have different lengths');
    end;
    Exit(MakeMatrixFromArray(rows));
  end
  else
  begin
    // parse as vector: comma-separated scalar expressions
    SetLength(rowVec, 0);
    while True do
    begin
      tmpVal := ParserExpression;
      if tmpVal.kind <> vkScalar then
        raise Exception.Create('Vector elements must be scalars');
      SetLength(rowVec, Length(rowVec)+1);
      rowVec[High(rowVec)] := tmpVal.scalar;
      SkipSpaces;
      if Peek = ',' then Next else Break;
    end;
    Expect(']');
    Exit(MakeVectorFromArray(rowVec));
  end;
end;

function TParserr.ParserFunction(const Name: string): TValue;
var arg1, arg2: TValue;
begin
  Expect('(');
  arg1 := ParserExpression;
  SkipSpaces;
  if Peek = ',' then
  begin
    Next;
    arg2 := ParserExpression;
  end;
  Expect(')');
  if SameText(Name, 'ln') then
  begin
    if arg1.kind <> vkScalar then raise Exception.Create('ln expects scalar');
    Result := MakeScalar(Ln(arg1.scalar));
    Exit;
  end;
  if SameText(Name, 'log') then
  begin
    if arg1.kind <> vkScalar then raise Exception.Create('log expects scalar');
    Result := MakeScalar(Log10(arg1.scalar));
    Exit;
  end;
  if SameText(Name, 'exp') then
  begin
    if arg1.kind <> vkScalar then raise Exception.Create('exp expects scalar');
    Result := MakeScalar(Exp(arg1.scalar));
    Exit;
  end;
  if SameText(Name, 'dot') then
  begin
    Result := DotProduct(arg1, arg2);
    Exit;
  end;
  if SameText(Name, 'matmul') then
  begin
    Result := MatMulFunction(arg1, arg2);
    Exit;
  end;
  raise Exception.Create('Unknown function: ' + Name);
end;

function TParserr.ParserPrimary: TValue;
begin
  SkipSpaces;
  case Peek of
    '(':
      begin
        Next; // (
        Result := ParserExpression;
        Expect(')');
      end;
    '[':
      begin
        Result := ParserArrayLiteral;
      end;
    '0'..'9', '.':
      Result := ParserNumber;
    '+', '-':
      begin
        // unary + or - handled in ParserUnary
        Result := ParserUnary;
      end;
  else
    if Peek in ['A'..'Z','a'..'z','_'] then Exit(ParserIdentifier);
    raise Exception.CreateFmt('Unexpected character "%s" at pos %d', [Peek, FPos]);
  end;
end;

function TParserr.ParserUnary: TValue;
var op: Char;
begin
  SkipSpaces;
  op := Peek;
  if (op = '+') or (op = '-') then
  begin
    Next;
    Result := ParserUnary;
    if op = '-' then
    begin
      case Result.kind of
        vkScalar: Result.scalar := -Result.scalar;
        vkVector: for var i := 0 to High(Result.vec) do Result.vec[i] := -Result.vec[i];
        vkMatrix: for var r := 0 to High(Result.mat) do for var c := 0 to High(Result.mat[r]) do Result.mat[r][c] := -Result.mat[r][c];
      end;
    end;
  end
  else
    Result := ParserPrimary;
end;

function TParserr.ParserFactor: TValue;
var base, exponent: TValue;
begin
  base := ParserUnary;
  SkipSpaces;
  while Peek = '^' do
  begin
    Next;
    exponent := ParserUnary;
    base := PowVal(base, exponent);
  end;
  Result := base;
end;

function TParserr.ParserTerm: TValue;
var left, right: TValue; op: Char;
begin
  left := ParserFactor;
  SkipSpaces;
  while Peek in ['*', '/', '%'] do
  begin
    op := Peek; Next;
    right := ParserFactor;
    case op of
      '*': left := MulVal(left, right);
      '/': left := DivVal(left, right);
      '%': begin
             // mod: only scalars or elementwise already supported via elementwise functions
             if (left.kind = vkScalar) and (right.kind = vkScalar) then left := MakeScalar(Fmod(left.scalar, right.scalar))
             else if (left.kind = vkVector) and (right.kind = vkVector) then left := MakeVectorFromArray(ElementwiseVecOp(left.vec, right.vec, '%'))
             else if (left.kind = vkMatrix) and (right.kind = vkMatrix) then left := MakeMatrixFromArray(ElementwiseMatOp(left.mat, right.mat, '%'))
             else raise Exception.Create('Unsupported % operation for types');
           end;
    end;
    SkipSpaces;
  end;
  Result := left;
end;

function TParserr.ParserAddSub: TValue;
var left, right: TValue; op: Char;
begin
  left := ParserTerm;
  SkipSpaces;
  while Peek in ['+', '-'] do
  begin
    op := Peek; Next;
    right := ParserTerm;
    case op of
      '+': left := AddVal(left, right);
      '-': left := SubVal(left, right);
    end;
    SkipSpaces;
  end;
  Result := left;
end;

function TParserr.ParserComparison: TValue;
var left, right: TValue; op: string;
begin
  left := ParserAddSub;
  SkipSpaces;
  // allow multi-char comparisons
  if Peek in ['<', '>', '='] then
  begin
    op := Peek; Next;
    if (op = '<') and (Peek in ['=', '>']) then begin op := op + Peek; Next; end;
    if (op = '>') and (Peek = '=') then begin op := op + Peek; Next; end;
    right := ParserAddSub;
    // comparisons operate on scalars; you can extend comparison semantics later
    if (left.kind <> vkScalar) or (right.kind <> vkScalar) then
      raise Exception.Create('Comparisons currently defined for scalars only');
    case op of
      '<': Result := MakeBool(left.scalar < right.scalar);
      '>': Result := MakeBool(left.scalar > right.scalar);
      '<=': Result := MakeBool(left.scalar <= right.scalar);
      '>=': Result := MakeBool(left.scalar >= right.scalar);
      '=': Result := MakeBool(Abs(left.scalar-right.scalar) < 1e-12);
      '<>': Result := MakeBool(Abs(left.scalar-right.scalar) >= 1e-12);
    else raise Exception.Create('Unknown comparison');
    end;
  end
  else Result := left;
end;

function TParserr.ParserLogic: TValue;
var left, right: TValue; word: string;
begin
  left := ParserComparison;
  SkipSpaces;
  while True do
  begin
    // parse tokens like 'and' 'or'
    word := '';
    while Peek in ['A'..'Z','a'..'z'] do begin word := word + Peek; Next; end;
    if word = '' then Break;
    word := LowerCase(word);
    SkipSpaces;
    right := ParserComparison;
    if (left.kind <> vkBoolean) or (right.kind <> vkBoolean) then
      raise Exception.Create('Logical operators expect boolean operands');
    if word = 'and' then left := MakeBool(left.bool and right.bool)
    else if word = 'or' then left := MakeBool(left.bool or right.bool)
    else raise Exception.Create('Unknown logical operator: ' + word);
    SkipSpaces;
  end;
  Result := left;
end;

function TParserr.ParserExpression: TValue;
begin
  Result := ParserLogic;
end;

function TParserr.ParserStatement: TValue;
var name: string; found: Boolean; rhs: TValue;
    // helper to rewind if not assignment
    procedure Rewind(n: Integer);
    begin
      FPos := Max(1, FPos - n);
    end;
var startPos: Integer;
begin
  SkipSpaces;
  if Peek in ['A'..'Z','a'..'z','_'] then
  begin
    // read an identifier but don't consume permanently until we confirm assignment
    startPos := FPos;
    name := '';
    while Peek in ['A'..'Z','a'..'z','0'..'9','_'] do
    begin
      name := name + Peek; Next;
    end;
    SkipSpaces;
    if (Peek = ':') or (Peek = '=') then
    begin
      // assignment
      if Peek = ':' then
      begin
        Next; Expect('=');
      end
      else Next; // consume '='
      rhs := ParserExpression;
      VarSet(name, rhs);
      Exit(rhs);
    end
    else
    begin
      // not assignment: rewind to start and parse expression normally
      FPos := startPos;
    end;
  end;
  Result := ParserExpression;
end;

function TParserr.ParserTop: TValue;
begin
  FPos := 1;
  Result := ParserStatement;
  SkipSpaces;
  if Peek <> #0 then
    raise Exception.CreateFmt('Unexpected trailing characters at pos %d', [FPos]);
end;

{ ---------- External entry points ---------- }

procedure EvalExpression(const Expr: string);
var parser: TParserr; value: TValue;
begin
  parser := TParserr.Create;
  try
    parser.FExpr := Expr;
    parser.FPos := 1;
    value := parser.ParserTop;
    case value.kind of
      vkScalar: Writeln('= ', value.scalar:0:8);
      vkBoolean: Writeln('= ', BoolToStr(value.bool, True));
      vkVector:
        begin
          Write('= [');
          for var i := 0 to High(value.vec) do
          begin
            if i>0 then Write(', ');
            Write(value.vec[i]:0:6);
          end;
          Writeln(']');
        end;
      vkMatrix:
        begin
          Writeln('= [');
          for var r := 0 to High(value.mat) do
          begin
            Write('  [');
            for var c := 0 to High(value.mat[r]) do
            begin
              if c>0 then Write(', ');
              Write(value.mat[r][c]:0:6);
            end;
            Writeln(']');
          end;
          Writeln(']');
        end;
    end;
  finally
    parser.Free;
  end;
end;

function EvalToValue(const Expr: string): TValue;
var parser: TParserr;
begin
  parser := TParserr.Create;
  try
    parser.FExpr := Expr;
    parser.FPos := 1;
    Result := parser.ParserTop;
  finally
    parser.Free;
  end;
end;

{ ---------- Initialization / finalization ---------- }

initialization
  Vars := TStringList.Create;

finalization
  Vars.Free;

end.

end.

