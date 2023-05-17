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
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

KMatrix<Rational> LatticeMinpoly(BigCplx alpha, BigCplx beta, int d, int O)
{
    var pi = BigReal.Pi(alpha.O);
    var N = new Rational(BigInteger.Pow(10, O));
    var mat = new KMatrix<Rational>(Rational.KZero(), d, d).Zero;
    var ai = alpha.One;
    for (int i = 0; i < mat.M - 1; i++)
    {
        mat.Coefs[i, i] = Rational.KOne();
        var aipi = ai.RealPart + pi * ai.ImaginaryPart; // Re(𝛼^i) + π * Im(𝛼^i)
        mat.Coefs[i, mat.N - 1] = (aipi.ToRational * N).RoundEven;
        ai *= alpha;
    }

    var bpi = beta.RealPart + pi * beta.ImaginaryPart; // Re(β) + π * Im(β)
    mat.Coefs[mat.N - 1, mat.N - 1] = (bpi.ToRational * N).RoundEven;
    return mat;
}

Rational[] Run_fpLLL(KMatrix<Rational> mat, int O1)
{
    var matStr = string.Format("[{0}]", mat.Rows.Select(r => r.Glue(" ")).Glue(fmt: "[{0}]"));
    File.WriteAllText("tmp", matStr);
    
    var process = new Process();
    process.StartInfo.FileName = "fplll";
    process.StartInfo.Arguments = $"-a lll -f mpfr -z mpz -m proved -p {O1} -of uk tmp";
    process.StartInfo.RedirectStandardInput = false;
    process.StartInfo.RedirectStandardOutput = true;
    process.StartInfo.UseShellExecute = false;
    process.Start();
    process.WaitForExit();
    var res = process.StandardOutput.ReadToEnd();
    File.Delete("tmp");
    process.Close();
    process.Dispose();
    
    var coefs = res.Replace("\n", "").Replace("[[", "").Replace("]]", "")
        .Split("],[")[0].Split(",").Select(e => new Rational(-BigInteger.Parse(e)))
        .ToArray();
    
    return coefs;
}

void fpLLLminPoly(int r, int s, int O)
{
    GlobalStopWatch.AddLap();
    var n = r * s; // Expected polynomial degree
    var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, s, O); // a = 3^(1/r) - 2^(1/s)
    var beta = alpha.Pow(n);
    var mat = LatticeMinpoly(BigCplx.FromBigReal(alpha), BigCplx.FromBigReal(beta), n + 1, O);

    var coefs = Run_fpLLL(mat, O).ToKMatrix();
    Console.WriteLine(coefs);
    GlobalStopWatch.Show($"fpLLL min poly a = 3^(1/{r}) - 2^(1/{s})");
    var P = FG.KPoly('X', coefs.ToArray()).Monic;
    Console.WriteLine($"P = {P}");
    Console.WriteLine("P(alpha) = 0 : {0}", P.Substitute(alpha).ToBigReal(O - 4).IsZero());
    Console.WriteLine();
}

void arbQAut(KPoly<Rational> P, int O1, int O2)
{
    var Pc = P.ToBcPoly(O2);
    GlobalStopWatch.AddLap();
    var cRoots = Run_arbPolyRoots(P, O2);
    var (alpha, beta) = cRoots.Grid2D().Where(e => !e.t1.Equals(e.t2) && !e.t1.ToBigCplx(O1).Equals(-e.t2.ToBigCplx(O1))).First();
    Console.WriteLine(cRoots.Aggregate(Pc.One, (acc, c) => acc * (Pc.X - c)).ToBcPoly(O1).ToIntPoly());
    Console.WriteLine(P);
    GlobalStopWatch.Show("arb Roots");
    GlobalStopWatch.AddLap();
    var mat = LatticeMinpoly(alpha, beta, P.Degree + 1, O1);
    GlobalStopWatch.Show("Lattice");
    GlobalStopWatch.AddLap();
    var coefs = Run_fpLLL(mat, O1);
    GlobalStopWatch.Show($"fpLLL");
    
    var P0 = new KPoly<Rational>('x', Rational.KZero(), coefs.SkipLast(1).ToArray());
    var fact = (P0.Substitute(alpha) / beta).ToBigCplx(O1);
    P0 /= fact.RealPart.ToRational.RoundEven;
    
    var y = FG.EPoly(P, 'y');
    var P1 = P0.Substitute(y);
    Console.WriteLine($"With P = {P} and P(y) = 0");
    Console.WriteLine($"KAut : y -> {P1}");
    var gr = Group.KAut(P1);
    DisplayGroup.HeadElements(gr);
    Console.WriteLine();
}

BigCplx[] ReadStr(string[] rootsStr, int O)
{
    var cRoots = new List<BigCplx>();
    foreach (var r in rootsStr)
    {
        if (string.IsNullOrEmpty(r))
            continue;
        
        if (r.Contains(" + ") || r.Contains(" - "))
        {
            var s0 = r.Split(" + ");
            var i = 1;
            if (r.Contains(" - "))
            {
                s0 = r.Split(" - ");
                i = -1;
            }

            var (re, im) = (s0[0], s0[1].Replace("*I", ""));
            var bre = BigReal.FromString(re, O);
            var bim = BigReal.FromString(im, O);
            var rc = BigCplx.FromBigReal(bre, i * bim);
            cRoots.Add(rc);
        }
        else
        {
            var s0 = r.Replace("*I", "");
            var b = BigReal.FromString(s0, O);
            var rc = !r.Contains("*I") ? BigCplx.FromBigReal(b) : BigCplx.FromBigReal(b.Zero, b);
            cRoots.Add(rc);
        }
    }

    return cRoots.ToArray();
}

BigCplx[] Run_arbPolyRoots(KPoly<Rational> P, int O)
{
    if (P.Coefs.Any(c => !c.Denom.IsOne))
        throw new();
    
    var process = new Process();
    // assuming the gist https://gist.github.com/aidevnn/c60cdf55d2fe84a45caeb0853ed03aef
    // is already compiled and the binary is accessible 
    
    process.StartInfo.FileName = "polyroots"; 
    process.StartInfo.Arguments = string.Format("-refine {0} coeffs {1}", O, P.Coefs.Glue(" "));
    process.StartInfo.RedirectStandardInput = false;
    process.StartInfo.RedirectStandardOutput = true;
    process.StartInfo.UseShellExecute = false;
    process.Start();
    process.WaitForExit();
    var res = process.StandardOutput.ReadToEnd().Split("\n");
    process.Close();
    process.Dispose();

    return ReadStr(res, O);
}

{
    GlobalStopWatch.Restart();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    
    // MinPoly of a = 3^(1/r) - 2^(1/s)
    fpLLLminPoly(r: 2, s: 2, O: 20);
    fpLLLminPoly(r: 3, s: 3, O: 40);
    fpLLLminPoly(r: 2, s: 5, O: 50);
    fpLLLminPoly(r: 3, s: 4, O: 70);
    fpLLLminPoly(r: 2, s: 7, O: 90);
    fpLLLminPoly(r: 3, s: 5, O: 90);
    fpLLLminPoly(r: 4, s: 4, O: 90);
    fpLLLminPoly(r: 4, s: 5, O: 120);
    fpLLLminPoly(r: 5, s: 5, O: 160);
}

{
    var x = FG.QPoly();
    arbQAut(x.Pow(6) + 108, O1: 20, O2: 30);
    arbQAut(x.Pow(10) + 10 * x.Pow(8) + 125 * x.Pow(6) + 500 * x.Pow(4) + 2500 * x.Pow(2) + 4000, O1: 45, O2: 60);
    arbQAut(x.Pow(12) + 96 * x.Pow(8) + 1664 * x.Pow(6) - 16128 * x.Pow(4) + 165888 * x.Pow(2) + 331776, O1: 110, O2: 130);
    arbQAut(x.Pow(18) + 171 * x.Pow(12) + 5130 * x.Pow(6) + 27, O1: 120, O2: 140);
    arbQAut(x.Pow(20) + 2500 * x.Pow(10) + 50000, O1: 100, O2: 120);
    
    var P = x.Pow(21) - 84 * x.Pow(19) + 2436 * x.Pow(17) - 31136 * x.Pow(15) + 2312 * x.Pow(14) + 203840 * x.Pow(13) -
        30688 * x.Pow(12) - 733824 * x.Pow(11) + 152992 * x.Pow(10) + 1480192 * x.Pow(9) - 359296 * x.Pow(8) -
        1628096 * x.Pow(7) + 413952 * x.Pow(6) + 892416 * x.Pow(5) - 225792 * x.Pow(4) - 189952 * x.Pow(3) + 50176 * x.Pow(2) +
        3584 * x - 512; // C7x:C3
    arbQAut(P, O1: 250, O2: 300);
}
