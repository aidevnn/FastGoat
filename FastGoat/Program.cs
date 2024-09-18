using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using System.Reflection.Emit;
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
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

BigReal GenBR((int M0, int M1) m, (int E0, int E1) e, int O)
{
    var a0 = Rng.Next(m.M0, m.M1 + 1);
    var e0 = Rng.Next(e.E0, e.E1 + 1);
    return BigReal.FromBigInteger(a0, O).Mul10PowN(e0);
}

string fmt(double a, int O) => String.Format($"{{0:E{O}}}", a);

void BigRealAdd()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range()
        .Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), GenBR((-M, M), (-E, E), O)))
        .Select(e => Rng.NextDouble() < 0.025 ? (e.Item1, e.Item1) : Rng.NextDouble() < 0.5 ? e : (e.Item2, e.Item1))
        .ToArray();
    GlobalStopWatch.Restart();
    var err = Double.Pow(10, -O + 2);
    foreach (var (a, b) in lt)
    {
        var a_b = a + b;
        var a_b2 = a.ToDouble + b.ToDouble;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if (diff > err)
        {
            Console.WriteLine(new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b, bd = $"{fmt(b.ToDouble, O)}", b0 = b.Details });
            Console.WriteLine(
                new
                {
                    a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealAdd");
}

void BigRealMul()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range()
        .Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), GenBR((-M, M), (-E, E), O)))
        .Select(e => Rng.NextDouble() < 0.025 ? (e.Item1, e.Item1) : Rng.NextDouble() < 0.5 ? e : (e.Item2, e.Item1))
        .ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var (a, b) in lt)
    {
        var a_b = a * b;
        var a_b2 = a.ToDouble * b.ToDouble;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if (diff > err)
        {
            Console.WriteLine(new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b, bd = $"{fmt(b.ToDouble, O)}", b0 = b.Details });
            Console.WriteLine(
                new
                {
                    a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealMul");
}

void BigRealKMul()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range().Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), Rng.Next(-M, M))).ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var (a, b) in lt)
    {
        var a_b = a * b;
        var a_b2 = a.ToDouble * b;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if (diff > err)
        {
            Console.WriteLine(new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b, b0 = BigReal.FromBigInteger(b, O) });
            Console.WriteLine(
                new
                {
                    a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealKMul");
}

void BigRealDiv()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range()
        .Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), GenBR((-M, M), (-E, E), O)))
        .Select(e => Rng.NextDouble() < 0.025 ? (e.Item1, e.Item1) : Rng.NextDouble() < 0.5 ? e : (e.Item2, e.Item1))
        .ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var (a, b) in lt)
    {
        if (b.IsZero())
            continue;
        
        var a_b = a / b;
        var a_b2 = a.ToDouble / b.ToDouble;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if (diff > err)
        {
            Console.WriteLine(new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b });
            Console.WriteLine(
                new
                {
                    a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealDiv");
}

void BigRealStr()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range().Select(_ => GenBR((-M * M, M * M), (-2 * E, 2 * E), O)).ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var a in lt)
    {
        var a_b = a;
        var a_b2 = Double.Parse($"{a}");
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if (diff > err)
        {
            Console.WriteLine(new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details });
            Console.WriteLine(
                new
                {
                    a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealStr");
}

void BigRealBInt()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range().Select(_ => Rng.Next(-M * M, M * M) * BigInteger.Pow(10, Rng.Next(0, E))).ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var a in lt)
    {
        var a_b = BigReal.FromBigInteger(a, O);
        var a_b2 = Double.Parse($"{a}");
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if (diff > err)
        {
            Console.WriteLine(
                new
                {
                    a, a_b.Details, a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealBInt");
}

void BigRealRat()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range().Select(_ => GenBR((-M * M, M * M), (-2 * E, 2 * E), O)).ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var a in lt)
    {
        var a_b = a;
        var a_b2 = (double)a.ToRational;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if (diff > err)
        {
            Console.WriteLine(
                new
                {
                    a, a_b.Details, a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealRat");
}

{
    GlobalStopWatch.Restart();
    GlobalStopWatch.AddLap();
    
    BigRealBInt();
    BigRealStr();
    BigRealAdd();
    BigRealMul();
    BigRealKMul();
    BigRealDiv();
    BigRealRat();
    
    GlobalStopWatch.Show("END");
}