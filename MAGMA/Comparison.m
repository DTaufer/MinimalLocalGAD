// load all functions first

R<x,y,z> := PolynomialRing(Rationals(),3);
CoordinateChange := [Random([-10..10])*x+Random([-10..10])*y+Random([-10..10])*z : i in [1..Rank(R)]];

F := x^2*y + x*y*z + y^3;
time for i in [1..100] do ms:=MinimalSupportsA(F,4); end for; // 16.05
time for i in [1..100] do ms:=MinimalSupportsB(F,4); end for; // 0.82
time for i in [1..100] do ms:=MinimalSupportsC(F,4); end for; // 0.03

time for i in [1..100] do ms:=MinimalSupportsA(Evaluate(F,CoordinateChange),4); end for; // 55.82
time for i in [1..100] do ms:=MinimalSupportsB(Evaluate(F,CoordinateChange),4); end for; // 7.15
time for i in [1..100] do ms:=MinimalSupportsC(Evaluate(F,CoordinateChange),4); end for; // 5.23


F := x^2*y*z;
time for i in [1..10] do ms:=MinimalSupportsA(F,4); end for; // 10.53
time for i in [1..10] do ms:=MinimalSupportsB(F,4); end for; // 0.84
time for i in [1..10] do ms:=MinimalSupportsC(F,4); end for; // 0.01

time for i in [1..10] do ms:=MinimalSupportsA(Evaluate(F,CoordinateChange),4); end for; // 14.13
time for i in [1..10] do ms:=MinimalSupportsB(Evaluate(F,CoordinateChange),4); end for; // 3.34
time for i in [1..10] do ms:=MinimalSupportsC(Evaluate(F,CoordinateChange),4); end for; // 3.18


F := x^2*(y+z) + y^2*(x+z) + z^2*(x+y);
time for i in [1..20] do ms:=MinimalSupportsA(F,5); end for; // 11.18
time for i in [1..20] do ms:=MinimalSupportsB(F,5); end for; // 0.50
time for i in [1..20] do ms:=MinimalSupportsC(F,5); end for; // 0.44

time for i in [1..20] do ms:=MinimalSupportsA(Evaluate(F,CoordinateChange),5); end for; // 60.50
time for i in [1..20] do ms:=MinimalSupportsB(Evaluate(F,CoordinateChange),5); end for; // 25.11
time for i in [1..20] do ms:=MinimalSupportsC(Evaluate(F,CoordinateChange),5); end for; // 13.74


F := x^2*y^2*z;
time for i in [1..5] do ms:=MinimalSupportsA(F,6); end for; // NO
time for i in [1..5] do ms:=MinimalSupportsB(F,6); end for; // 17.98
time for i in [1..5] do ms:=MinimalSupportsC(F,6); end for; // 0.01

time for i in [1..5] do ms:=MinimalSupportsA(Evaluate(F,CoordinateChange),6); end for; // NO
time for i in [1..5] do ms:=MinimalSupportsB(Evaluate(F,CoordinateChange),6); end for; // 113.35
time for i in [1..5] do ms:=MinimalSupportsC(Evaluate(F,CoordinateChange),6); end for; // 54.64


R<x,y,z,u> := PolynomialRing(Rationals(),4);
CoordinateChange := [Random([-5..5])*x+Random([-5..5])*y+Random([-5..5])*z+Random([-5..5])*u : i in [1..Rank(R)]];

F := x^3 + y^3 + x^2*u + x*u*z + u*z^2;
time for i in [1..3] do ms:=MinimalSupportsA(F,7); end for; // NO
time for i in [1..3] do ms:=MinimalSupportsB(F,7); end for; // 21.34
time for i in [1..3] do ms:=MinimalSupportsC(F,7); end for; // 0.97

time for i in [1..3] do ms:=MinimalSupportsA(Evaluate(F,CoordinateChange),7); end for; // NO
time for i in [1..3] do ms:=MinimalSupportsB(Evaluate(F,CoordinateChange),7); end for; // NO
time for i in [1..3] do ms:=MinimalSupportsC(Evaluate(F,CoordinateChange),7); end for; // NO


F := x*(x^3 + x^2*y + x*z^2 + y^3 + z^3 + u^3);
time for i in [1..3] do ms:=MinimalSupportsA(F,8); end for; // NO
time for i in [1..3] do ms:=MinimalSupportsB(F,8); end for; // 143.63
time for i in [1..3] do ms:=MinimalSupportsC(F,8); end for; // 21.15

time for i in [1..3] do ms:=MinimalSupportsA(Evaluate(F,CoordinateChange),8); end for; // NO
time for i in [1..3] do ms:=MinimalSupportsB(Evaluate(F,CoordinateChange),8); end for; // NO
time for i in [1..3] do ms:=MinimalSupportsC(Evaluate(F,CoordinateChange),8); end for; // NO
