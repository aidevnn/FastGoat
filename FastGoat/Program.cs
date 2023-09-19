using System.Diagnostics;
using System.IO.IsolatedStorage;
using System.Numerics;
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

// H.E.Rose, A Course on FInite Groups, Problem 3.22 page 64
ConcreteGroup<KMatrix<Cnf>> DiCyclic(int m, bool validation = false)
{
    if (m < 2)
        throw new GroupException(GroupExceptionType.BaseGroup);

    var w = Cnf.Nth(2 * m);
    var gl = FG.GLnK("C", 2, w);
    var Am = gl[w, 0, 0, w.Inv()];
    var B = gl[0, 1, -1, 0];
    var I2 = gl.Neutral();
    if (validation && (!Am.Pow(m).Equals(-I2) || !B.Pow(2).Equals(-I2)))
        throw new();

    var dicm = Group.Generate($"Dic{m}gl", gl, Am, B);
    if (validation && !dicm.IsIsomorphicTo(FG.DiCyclic(m)))
        throw new();

    return dicm;
}

ConcreteGroup<KMatrix<EPoly<ZnInt>>> DiCyclicFp(int m, bool validation = false)
{
    if (m < 2)
        throw new GroupException(GroupExceptionType.BaseGroup);

    var p = Primes10000.First(p => Gcd(p, 2 * m) == 1);
    var w = FG.NthRootFp(2 * m, p).PrimitivesRoots().First();
    var gl = FG.GLnK($"F{p}[{w.Poly.x}/({w.F})]", 2, w);
    var Am = gl[w, 0, 0, w.Inv()];
    var B = gl[0, 1, -1, 0];
    var C = gl[-1, 0, 0, -1];
    if (validation && (!gl.Times(Am, m).Equals(C) || !gl.Times(B, 2).Equals(C)))
        throw new();

    var dicm = Group.Generate($"Dic{m}gl", gl, Am, B);
    if (validation && !dicm.IsIsomorphicTo(FG.DiCyclic(m)))
        throw new();

    return dicm;
}

(int p, int a, int ai) QuaternionGeneratorsSearch(int m){

    var (p, a) = Primes10000.Select(p => (p, Solve_k_pow_m_equal_one_mod_n_strict(p, 2 * m))).First(e => e.Item2 != -1);
    var ai = UnInvertible(p).First(k => k.Key == a).Value;
    return (p, a, ai);
}

ConcreteGroup<Mat> DiCyclicZnInt(int m, bool validation = false)
{
    if (m < 2)
        throw new GroupException(GroupExceptionType.BaseGroup);

    var (p, a, ai) = QuaternionGeneratorsSearch(m);
    var gl = new GL(2, p);
    var Am = gl[a, 0, 0, ai];
    var B = gl[0, 1, -1, 0];
    var C = gl[-1, 0, 0, -1];
    if (validation && (!gl.Times(Am, m).Equals(C) || !gl.Times(B, 2).Equals(C)))
        throw new();

    var dicm = Group.Generate($"Dic{m}gl", gl, Am, B);
    if (validation && !dicm.IsIsomorphicTo(FG.DiCyclic(m)))
        throw new();

    return dicm;
}

void BenchDic()
{
    GlobalStopWatch.Restart();
    for (int m = 2; m < 30; m++)
    {
        for (int i = 0; i < 3; i++)
        {
            GlobalStopWatch.AddLap();
            GlobalStopWatch.Show($"{DiCyclicZnInt(m)}");
            GlobalStopWatch.AddLap();
            GlobalStopWatch.Show($"{FG.DiCyclicSdp(m)}");
            Console.WriteLine();
        }

        Console.WriteLine();
    }
}

void BenchQuaternion()
{
    GlobalStopWatch.Restart();
    for (int k = 1; k <= 5; k++)
    {
        var m = 2.Pow(k);
        for (int i = 0; i < 3; i++)
        {
            GlobalStopWatch.AddLap();
            GlobalStopWatch.Show($"{DiCyclicZnInt(m).ShortName}");
            GlobalStopWatch.AddLap();
            GlobalStopWatch.Show($"{FG.DiCyclicSdp(m).ShortName}");
            Console.WriteLine();
        }

        Console.WriteLine();
    }
}

void BenchQuaternionGensSearch()
{
    for (int k = 1; k <= 7; k++)
    {
        var m = 2.Pow(k);
        for (int i = 0; i < 3; i++)
        {
            GlobalStopWatch.AddLap();
            var s = QuaternionGeneratorsSearch(m);
            GlobalStopWatch.Show($"Search Q{m * 4} => p={s.p} a={s.a} ai={s.ai}");
        }
    }
}

void DiCyclicGLnC()
{
    for (int m = 2; m < 11; ++m)
    {
        var Qm = DiCyclic(m, validation: false);
        DisplayGroup.HeadOrders(Qm);
        foreach (var mat in Qm.GetGenerators())
        {
            Console.WriteLine($"Gen order {Qm.ElementsOrders[mat]}");
            Console.WriteLine(mat);
        }
        
        Console.WriteLine();
    }
}

void DiCyclicGLnq(bool validation = false)
{
    for (int m = 2; m < 21; ++m)
    {
        var Qm = DiCyclicFp(m, validation);
        DisplayGroup.HeadOrders(Qm);
        foreach (var mat in Qm.GetGenerators())
        {
            Console.WriteLine($"Gen order {Qm.ElementsOrders[mat]}");
            Console.WriteLine(mat);
        }

        Console.WriteLine();
    }
}

void DiCyclicGLnp(bool validation = false)
{
    for (int m = 2; m < 31; ++m)
    {
        var Qm = DiCyclicZnInt(m, validation);
        DisplayGroup.HeadOrders(Qm);
        foreach (var mat in Qm.GetGenerators())
        {
            Console.WriteLine($"Gen order {Qm.ElementsOrders[mat]}");
            Console.WriteLine(mat);
        }

        Console.WriteLine();
    }
}

void QuaternionGroupGeneratorsZpInt()
{
    for (int n = 1; n < 8; ++n)
    {
        var m = 2.Pow(n + 2) / 4;
        var Qm = DiCyclicZnInt(m);
        Qm.SetName($"Q{m * 4}");
        DisplayGroup.HeadOrders(Qm);
        foreach (var mat in Qm.GetGenerators())
        {
            Console.WriteLine($"Gen order {Qm.ElementsOrders[mat]}");
            Console.WriteLine(mat);
        }
        
        Console.WriteLine();
    }
}

{
    // BenchDic();
    // BenchQuaternion();
    // BenchQuaternionGensSearch();

    DiCyclicGLnC();
    DiCyclicGLnq();
    DiCyclicGLnp();
    QuaternionGroupGeneratorsZpInt();
}

/* 
   |Q128| = 128
   Type        NonAbelianGroup
   BaseGroup   GL(2,193)
   
   Elements Orders : [1]:1, [2]:1, [4]:66, [8]:4, [16]:8, [32]:16, [64]:32
   
   Gen order 64
   [182,   0]
   [  0,  35]
   Gen order 4
   [  0,   1]
   [192,   0]
   
   |Q256| = 256
   Type        NonAbelianGroup
   BaseGroup   GL(2,257)
   
   Elements Orders : [1]:1, [2]:1, [4]:130, [8]:4, [16]:8, [32]:16, [64]:32, [128]:64
   
   Gen order 128
   [248,   0]
   [  0,  57]
   Gen order 4
   [  0,   1]
   [256,   0]
   
   |Q512| = 512
   Type        NonAbelianGroup
   BaseGroup   GL(2,257)
   
   Elements Orders : [1]:1, [2]:1, [4]:258, [8]:4, [16]:8, [32]:16, [64]:32, [128]:64, [256]:128
   
   Gen order 256
   [254,   0]
   [  0, 171]
   Gen order 4
   [  0,   1]
   [256,   0]
*/