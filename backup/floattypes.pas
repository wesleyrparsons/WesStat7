unit FloatTypes;

  {$mode ObjFPC}{$H+}

interface

// === Precision Selection ===
// Uncomment ONE of the following or define via compiler command line:
// {$DEFINE USE_SINGLE}
  {$DEFINE USE_DOUBLE}
// {$DEFINE USE_EXTENDED}

type
  {$IFDEF USE_SINGLE}
    TFloat = Single;
    const Epsilon: TFloat = 1e-6;
  {$ENDIF}

  {$IFDEF USE_DOUBLE}
    TFloat = Double;
    const Epsilon: TFloat = 1e-12;
  {$ENDIF}

  {$IFDEF USE_EXTENDED}
    TFloat = Extended;
    const Epsilon: TFloat = 1e-18;
  {$ENDIF}

type
      TVector = array of TFloat;
      TMatrix = array of TVector;
      TIVector = array of Integer;
      TIMatrix = array of TIVector;
      TFloat = Double;
      MessageType = array[0..50] of String[80];

implementation

end.
