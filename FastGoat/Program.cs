using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
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

{
    var d16 = FG.Dihedral(8);
    var d16a = d16.ToTable().gt;
    DisplayGroup.HeadElements(d16a);
    var subgroups = new AllSubgroups<Perm>(d16);
    var subgroupsa = new AllSubgroups<TableElt>(d16a);
    var subgroupsb = subgroups.ToTable();
    Console.WriteLine(subgroups.Infos);
    Console.WriteLine(subgroupsa.Infos);
    Console.WriteLine(subgroupsb.Infos);
}
