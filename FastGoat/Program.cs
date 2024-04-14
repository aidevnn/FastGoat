using System.Collections;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

WordGroup[] GetGroups(int ord)
{
    return GroupExt.DB.Select(s => s.Split(';'))
        .Where(s => int.Parse(s[0]) == ord)
        .Select(s =>
        {
            Logger.Level = LogLevel.Off;
            var g = FG.WordGroup(s[1], s[2]);
            Logger.Level = LogLevel.Level1;
            return g;
        })
        .ToArray();
}

void MatrixFormGroupsOrd8()
{
    var gl2 = FG.GLnq(2, 5);
    var so3 = FG.GO3p(3);
    foreach (var g8 in GetGroups(8))
    {
        Console.WriteLine(g8.ShortName);
        if (g8.Name == "Q8")
        {
            var iso = FG.Quaternion(8);
            DisplayGroup.Generators(iso);
            DisplayGroup.Generators(Group.IsomorphicSubgroup(gl2, g8));
        }
        else
        {
            if (g8.Name != "C2 x C2 x C2")
                DisplayGroup.Generators(Group.IsomorphicSubgroup(gl2, g8));

            if (g8.Name != "C8")
                DisplayGroup.Generators(Group.IsomorphicSubgroup(so3, g8));
        }
    }
}

void MatrixFormGroupsOrd9()
{
    var gl2 = new GL(2, 19);
    var a = FG.FqX(19);
    var e3 = a.Pow(6)[0].K;
    var e9 = a.Pow(2)[0].K;
    var gen9 = gl2[e9, 0, 0, 1];
    var gen3a = gl2[e3, 0, 0, 1];
    var gen3b = gl2[1, 0, 0, e3];
    var g9 = Group.Generate("C9mat", gl2, gen9);
    var g33 = Group.Generate("(C3 x C3)mat", gl2, gen3a, gen3b);
    var ab9 = FG.Abelian(9);
    var ab33 = FG.Abelian(3, 3);
    
    DisplayGroup.Generators(g9);
    DisplayGroup.AreIsomorphics(g9, ab9);
    Console.WriteLine();
    
    DisplayGroup.Generators(g33);
    DisplayGroup.AreIsomorphics(g33, ab33);
    Console.WriteLine();
}

void MatrixFormGroupsOrd10()
{
    var gl2 = FG.GLnp(2, 11);
    foreach (var g10 in GetGroups(10))
    {
        var iso = Group.IsomorphicSubgroup(gl2, g10);
        DisplayGroup.Generators(iso);
    }
}

void MatrixFormGroupsOrd12()
{
    var gl2 = FG.GLnp(2, 7);
    DisplayGroup.HeadOrders(gl2);
    foreach (var g12 in GetGroups(12))
    {
        Console.WriteLine(g12.ShortName);
        if (g12.Name == "A4")
        {
            var so33 = FG.SO3q(3);
            var iso = Group.IsomorphicSubgroup(so33, g12);
            DisplayGroup.Generators(iso);
        }
        else
        {
            var iso = Group.IsomorphicSubgroup(gl2, g12);
            DisplayGroup.Generators(iso);
        }
    }
}

void MatrixFormGroupsOrd14()
{
    var gl2 = FG.GLnq(2, 8);
    foreach (var g14 in GetGroups(14))
    {
        var iso = Group.IsomorphicSubgroup(gl2, g14);
        DisplayGroup.Generators(iso);
    }
}

{
    MatrixFormGroupsOrd8();
    MatrixFormGroupsOrd9();
    MatrixFormGroupsOrd10();
    MatrixFormGroupsOrd12();
    MatrixFormGroupsOrd14();
}