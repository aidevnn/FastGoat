
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void SU3_2()
{
    var gl34 = new GLnq(3, 4);
    var x = gl34.Fq['x'];

    var a = gl34[
        1, x, x,
        0, 1, x + 1,
        0, 0, 1
    ];

    var b = gl34[
        x, 1, 1,
        1, 1, 0,
        1, 0, 0
    ];

    var su3_2 = Group.Generate("SU3(2)", gl34, a, b);
    DisplayGroup.Head(su3_2);

    try
    {
        su3_2.AllSubgroups().DisplayAllSeries(); // hash collision
    }
    catch (Exception e)
    {
        Console.WriteLine(e);
    }
    
    var (su3_2p, gens_su3_2) = su3_2.ToPermGroup();
    DisplayGroup.Head(su3_2p);
    var allSubGr_su3_2p = su3_2p.AllSubgroups();
    allSubGr_su3_2p.DisplayAllSeries();
    
    var z = Group.Zentrum(su3_2);
    DisplayGroup.Head(z);
    var u3_2 = su3_2.Over(z, "U3(2)");
    DisplayGroup.Head(u3_2);
    var (u3_2p, gens_u3_2) = u3_2.ToPermGroup();
    DisplayGroup.Head(u3_2p);
    var allSubGr_u3_2p = u3_2p.AllSubgroups();
    allSubGr_u3_2p.DisplayAllSeries();
}

{
    var d8 = FG.DihedralWg(4);
    DisplayGroup.HeadElements(d8);
    var (d8p, mapGens) = d8.ToPermGroup();
    DisplayGroup.HeadElements(d8p);
    mapGens.Println();

    SU3_2();
}