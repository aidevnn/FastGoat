using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
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

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

// Pairs of Generators for Matrix Groups. I
// arXiv:2201.09155v1 [math.GR] 23 Jan 2022
// http://arxiv.org/abs/2201.09155v1
// D. E. Taylor
// Department of Pure Mathematics
// The University of Sydney
// Australia 2006

void test1()
{
    var primes50 = Primes10000.Where(p => p < 50).ToArray();
    for (int n = 2; n <= 5; n++)
    {
        foreach (var p in primes50.Where(p => FG.GLnqOrder(n, p) < FG.MatrixGroupMaxOrder))
        {
            var glOrd = FG.GLnqOrder(n, p);
            var gl = FG.GLnp(n, p);
            var (a, b) = gl.GetGenerators().Deconstruct();
            Console.WriteLine(gl.ShortName);
            Console.WriteLine(a);
            Console.WriteLine();
            Console.WriteLine(b);
            Console.WriteLine();
            Console.WriteLine();
            if (gl.Count() != glOrd)
                throw new();
        }

        Console.WriteLine("####################################################################");

        foreach (var p in Primes10000.Where(p => FG.SLnqOrder(n, p) < FG.MatrixGroupMaxOrder))
        {
            var slOrd = FG.SLnqOrder(n, p);
            var sl = FG.SLnp(n, p);
            var (a, b) = sl.GetGenerators().Deconstruct();
            Console.WriteLine(sl.ShortName);
            Console.WriteLine(a);
            Console.WriteLine();
            Console.WriteLine(b);
            Console.WriteLine();
            if (sl.Count() != slOrd)
                throw new();
        }
        
        Console.WriteLine("####################################################################");
    }
}

void GLnq(int n, int q)
{
    GlobalStopWatch.AddLap();
    var glOrd = FG.GLnqOrder(n, q);
    var gl = FG.GLnq(n, q);
    var (a, b) = gl.GetGenerators().Deconstruct();
    Console.WriteLine(gl.ShortName);
    Console.WriteLine(a);
    Console.WriteLine();
    Console.WriteLine(b);
    Console.WriteLine();
    if (gl.Count() != glOrd)
        throw new();

    GlobalStopWatch.Show();
    Console.WriteLine();
}

void SLnq(int n, int q)
{
    GlobalStopWatch.AddLap();
    var slOrd = FG.SLnqOrder(n, q);
    var sl = FG.SLnq(n, q);
    var (a, b) = sl.GetGenerators().Deconstruct();
    Console.WriteLine(sl.ShortName);
    Console.WriteLine(a);
    Console.WriteLine();
    Console.WriteLine(b);
    Console.WriteLine();
    if (sl.Count() != slOrd)
        throw new();

    GlobalStopWatch.Show();
    Console.WriteLine();
}

void test2()
{
    GlobalStopWatch.Restart();
    foreach (var q in new[] { 4, 8, 9, 16 })
        GLnq(2, q);
    
    GLnq(3, 4);
        
    Console.WriteLine("####################################################################");

    foreach (var q in new[] { 4, 8, 9, 16, 25, 27, 32 })
        SLnq(2, q);
    
    SLnq(3, 4);
        
    Console.WriteLine("####################################################################");
}

{
    test1();
    // test2(); // slow TODO
}