using System.Diagnostics.Tracing;
using System.Numerics;
using System.Security.Principal;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void Dihedral8()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly('X');
    var P = x.Pow(4) + 2;
    var gal = GaloisApplicationsPart2.GaloisGroupChebotarev(P, detail: true);

    var minPol = IntFactorisation.SplittingField(P, true)[0].F.SubstituteChar('X');
    Console.WriteLine(minPol);

    var (X1, y1) = FG.EPolyXc(minPol, 'y');
    var def = IntFactorisation.DeflateByN(minPol, 2);

    var facts = IntFactorisation.AlgebraicFactors(def.nf.Substitute(X1), true);
    var roots = facts.SelectMany(f => IntFactorisation.AlgebraicRoots(f.Substitute(X1.Pow(2)), true)).ToList();

    var subFields = GaloisTheory.SubFields(roots, 3).ToArray();
    var extTowers = GaloisApplications.ExtensionsTower(subFields);
    GaloisApplications.GaloisCorrespondence(extTowers);
    DisplayGroup.AreIsomorphics(gal, subFields.OrderBy(gc => gc.SubGr.Count()).Last().SubGr.SuperGroup!);
}

void Alternate4()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly('X');
    var P = x.Pow(4) + 8 * x + 12;
    var gal = GaloisApplicationsPart2.GaloisGroupChebotarev(P, detail: true);

    var minPol = IntFactorisation.SplittingField(P, true)[0].F.SubstituteChar('X');
    Console.WriteLine(minPol);

    var (X1, y1) = FG.EPolyXc(minPol, 'y');
    var def = IntFactorisation.DeflateByN(minPol, 2);

    var facts = IntFactorisation.AlgebraicFactors(def.nf.Substitute(X1), true);
    var roots = facts.SelectMany(f => IntFactorisation.AlgebraicRoots(f.Substitute(X1.Pow(2)), true)).ToList();

    var subFields = GaloisTheory.SubFields(roots, 3).ToArray();
    var extTowers = GaloisApplications.ExtensionsTower(subFields);
    GaloisApplications.GaloisCorrespondence(extTowers);
    DisplayGroup.AreIsomorphics(gal, subFields.OrderBy(gc => gc.SubGr.Count()).Last().SubGr.SuperGroup!);
}

void Dihedral10()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly('X');
    var P = x.Pow(5) - 5 * x + 12;
    var gal = GaloisApplicationsPart2.GaloisGroupChebotarev(P, detail: true);

    var minPol = IntFactorisation.SplittingField(P, true)[0].F.SubstituteChar('X');
    Console.WriteLine(minPol);

    var (X1, y1) = FG.EPolyXc(minPol, 'y');
    var def = IntFactorisation.DeflateByN(minPol, 2);

    var facts = IntFactorisation.AlgebraicFactors(def.nf.Substitute(X1), true);
    var roots = facts.SelectMany(f => IntFactorisation.AlgebraicRoots(f.Substitute(X1.Pow(2)), true)).ToList();

    var subFields = GaloisTheory.SubFields(roots, 3).ToArray();
    var extTowers = GaloisApplications.ExtensionsTower(subFields);
    GaloisApplications.GaloisCorrespondence(extTowers);
    DisplayGroup.AreIsomorphics(gal, subFields.OrderBy(gc => gc.SubGr.Count()).Last().SubGr.SuperGroup!);
}

void Meta20()
{
    Ring.DisplayPolynomial = MonomDisplay.StarPowFct;
    var x = FG.QPoly('X');
    var P = x.Pow(5) + 2;
    var gal = GaloisApplicationsPart2.GaloisGroupChebotarev(P, detail: true);

    var minPol = IntFactorisation.SplittingField(P, details: true, onlyPositifs: true)[0].F.SubstituteChar('X');
    Console.WriteLine(minPol);

    var (X1, y1) = FG.EPolyXc(minPol, 'y');
    var def = IntFactorisation.DeflateByN(minPol, 5);

    var facts = IntFactorisation.AlgebraicFactors(def.nf.Substitute(X1), true);
    var roots = facts.SelectMany(f => IntFactorisation.AlgebraicRoots(f.Substitute(X1.Pow(5)), true)).ToList();

    var subFields = GaloisTheory.SubFields(roots, 3).ToArray();
    var extTowers = GaloisApplications.ExtensionsTower(subFields);
    GaloisApplications.GaloisCorrespondence(extTowers);
    DisplayGroup.AreIsomorphics(gal, subFields.OrderBy(gc => gc.SubGr.Count()).Last().SubGr.SuperGroup!);
    DisplayGroup.AreIsomorphics(gal, FG.Frobenius(20).First());
}

void Dihedral12()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly('X');
    var P = x.Pow(6) + 2;
    var gal = GaloisApplicationsPart2.GaloisGroupChebotarev(P, detail: true);

    var minPol = IntFactorisation.SplittingField(P, true)[0].F.SubstituteChar('X');
    Console.WriteLine(minPol);

    var (X1, y1) = FG.EPolyXc(minPol, 'y');
    var def = IntFactorisation.DeflateByN(minPol, 3);

    var facts = IntFactorisation.AlgebraicFactors(def.nf.Substitute(X1), true);
    var roots = facts.SelectMany(f => IntFactorisation.AlgebraicRoots(f.Substitute(X1.Pow(3)), true)).ToList();

    var subFields = GaloisTheory.SubFields(roots, 3).ToArray();
    var extTowers = GaloisApplications.ExtensionsTower(subFields);
    GaloisApplications.GaloisCorrespondence(extTowers);
    DisplayGroup.AreIsomorphics(gal, subFields.OrderBy(gc => gc.SubGr.Count()).Last().SubGr.SuperGroup!);
    DisplayGroup.AreIsomorphics(gal, FG.Dihedral(6));
}

void C6xC2()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly('X');
    var P = x.Pow(12) - x.Pow(6) + 1;
    // var gal = GaloisApplicationsPart2.GaloisGroupChebotarev(P, detail: true);

    var (X1, y1) = FG.EPolyXc(P, 'y');
    var def = IntFactorisation.DeflateByN(P, 3);

    var facts = IntFactorisation.AlgebraicFactors(def.nf.Substitute(X1), true);
    var roots = facts.SelectMany(f => IntFactorisation.AlgebraicRoots(f.Substitute(X1.Pow(3)), true)).ToList();

    var subFields = GaloisTheory.SubFields(roots, 3).ToArray();
    var extTowers = GaloisApplications.ExtensionsTower(subFields);
    GaloisApplications.GaloisCorrespondence(extTowers);
    DisplayGroup.AreIsomorphics(FG.Abelian(6, 2), subFields.OrderBy(gc => gc.SubGr.Count()).Last().SubGr.SuperGroup!);
}

void F18()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly('X');
    var P = x.Pow(6) + 3 * x.Pow(3) + 3;
    var gal = GaloisApplicationsPart2.GaloisGroupChebotarev(P, detail: true);

    var minPol = IntFactorisation.SplittingField(P, details: true)[0].F.SubstituteChar('X');
    Console.WriteLine(minPol);

    var (X1, y1) = FG.EPolyXc(minPol, 'y');
    var def = IntFactorisation.DeflateByN(minPol, 3);

    var facts = IntFactorisation.AlgebraicFactors(def.nf.Substitute(X1), true);
    var roots = facts.SelectMany(f => IntFactorisation.AlgebraicRoots(f.Substitute(X1.Pow(3)), true)).ToList();

    var subFields = GaloisTheory.SubFields(roots, 3).ToArray();
    var extTowers = GaloisApplications.ExtensionsTower(subFields);
    GaloisApplications.GaloisCorrespondence(extTowers);
    DisplayGroup.AreIsomorphics(gal, subFields.OrderBy(gc => gc.SubGr.Count()).Last().SubGr.SuperGroup!);
    DisplayGroup.AreIsomorphics(gal, Group.SemiDirectProd(FG.Abelian(3, 3), FG.Abelian(2)));
}

{
    Alternate4();
    F18();
    Dihedral8();
    Dihedral10();
    Dihedral12();
    C6xC2();
    // Meta20();
}