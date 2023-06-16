using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public static class ExternLibs
{
    public static KMatrix<Rational> Run_fpLLL(KMatrix<Rational> mat)
    {
        var matStr = string.Format("[{0}]", mat.Rows.Select(r => r.Glue(" ")).Glue(fmt: "[{0}]"));
        if (File.Exists("tmp"))
            File.Delete("tmp");
        
        File.WriteAllText("tmp", matStr);
        if (matStr.Length > 30000)
            throw new($"{(mat.Dim, matStr.Length)}");

        var process = new Process();
        process.StartInfo.FileName = "fplll";
        process.StartInfo.Arguments = $"-a lll -z mpz -m proved tmp";
        process.StartInfo.RedirectStandardInput = false;
        process.StartInfo.RedirectStandardOutput = true;
        process.StartInfo.UseShellExecute = false;
        process.Start();
        process.WaitForExit();
        var res = process.StandardOutput.ReadToEnd();
        File.Delete("tmp");
        process.Close();
        process.Dispose();

        Console.WriteLine("Done");
        var coefs = res.Replace("\n", "").Replace("[[", "").Replace("]]", "").Split("][").ToArray();
        var lll = coefs.SelectMany(l => l.Trim().Split(" ").Select(e => new Rational(BigInteger.Parse(e)))).ToKMatrix(coefs.Length);
        return lll.T;
    }

    public static KPoly<Rational> Run_polysGcd(KPoly<Rational> f, KPoly<Rational> g)
    {
        
        var process = new Process();
        process.StartInfo.FileName = "polysgcd";
        process.StartInfo.Arguments = $"f {f.Coefs.Glue(" ")} g {g.Coefs.Glue(" ")}";
        process.StartInfo.RedirectStandardInput = false;
        process.StartInfo.RedirectStandardOutput = true;
        process.StartInfo.UseShellExecute = false;
        process.Start();
        process.WaitForExit();
        var res = process.StandardOutput.ReadToEnd();
        var (s0, s1) = res.Split("  ").Deconstruct();
        
        var x = FG.QPoly();
        var P0 = s1.Split(" ").Select((c, i) => x.Pow(i) * new Rational(BigInteger.Parse(c))).Aggregate(x.Zero, (acc, c) => acc + c);
        process.Close();
        process.Dispose();
        return P0;
    }

}