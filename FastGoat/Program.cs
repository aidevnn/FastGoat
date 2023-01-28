using System.Globalization;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using System.Numerics;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.Arm;
using System.Security.Cryptography.X509Certificates;
using System.Xml;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    
    FG.CharactersTable(FG.Abelian(2)).DisplayCells();
    FG.CharactersTable(FG.Abelian(3)).DisplayCells();
    FG.CharactersTable(FG.Abelian(2, 2)).DisplayCells();
    FG.CharactersTable(FG.Abelian(2, 3)).DisplayCells();
    FG.CharactersTable(FG.Abelian(6)).DisplayCells();
    FG.CharactersTable(FG.Abelian(4)).DisplayCells();
    FG.CharactersTable(FG.Abelian(2, 4)).DisplayCells();
    FG.CharactersTable(FG.Abelian(2, 2, 2)).DisplayCells();
    FG.CharactersTable(FG.Abelian(2, 2, 3)).DisplayCells();

    FG.CharactersTable(FG.Dihedral(3)).DisplayCells();
    FG.CharactersTable(FG.Dihedral(4)).DisplayCells();
    FG.CharactersTable(FG.Quaternion(8)).DisplayCells();
    FG.CharactersTable(FG.Alternate(4)).DisplayCells();
    FG.CharactersTable(FG.Dihedral(5)).DisplayCells();
    FG.CharactersTable(FG.Dihedral(6)).DisplayCells();
    FG.CharactersTable(FG.DiCyclic(3)).DisplayCells();
    FG.CharactersTable(Group.SemiDirectProd(new Cn(7), new Cn(3))).DisplayCells();
    FG.CharactersTable(Group.SemiDirectProd(new Cn(4), new Cn(4))).DisplayCells();
    FG.CharactersTable(Group.SemiDirectProd(FG.Abelian(3, 3), new Cn(4))).DisplayCells();
}

/* TODO filling Character Table with solutions
 
 
|(C3 x C3) x: C4| = 36
Type        NonAbelianGroup
BaseGroup   C3 x C3 x C4

[1   1   1   1   1   1]
[1  -1   1   1   I  -I]
[4   0  x0  x1   0   0]
[4   0  x2  x3   0   0]
[1   1   1   1  -1  -1]
[1  -1   1   1  -I   I]
|(C3 x C3) x: C4| = 36
System
    4·x0 + 4·x2 + 4
    4·x1 + 4·x3 + 4
    x0x1 + x2x3 + 4
    x0² + x2² + -5
    x1² + x3² + -5
Solve
    x3² + x3 + -2
    x2 + x3 + 1
    x1 + x3 + 1
    x0 + -x3
|(C3 x C3) x: C4| = 36
[Class     1  2 3a 3b 4a 4b]
[ Size     1  9  4  4  9  9]
[                          ]
[  X.1     1  1  1  1  1  1]
[  X.2     1 -1  1  1  I -I]
[  X.3     4  0  #  #  0  0]
[  X.4     4  0  #  #  0  0]
[  X.5     1  1  1  1 -1 -1]
[  X.6     1 -1  1  1 -I  I]
*/