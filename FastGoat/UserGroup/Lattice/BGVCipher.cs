using System.Numerics;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

public readonly struct BGVCipher
{
    public KPoly<Rational> A { get; }
    public KPoly<Rational> B { get; }

    public BGVCipher(KPoly<Rational> a, KPoly<Rational> b)
    {
        (A, B) = (a, b);
    }
    
    public BGVCipher CoefsMod(Rational q) => (A.CoefsMod(q), B.CoefsMod(q));

    public void Show(string name = "")
    {
        if(!string.IsNullOrEmpty(name))
            Console.WriteLine(name);
        Console.WriteLine($"A:{A}");
        Console.WriteLine($"B:{B}");
    }

    public static implicit operator BGVCipher((KPoly<Rational> a, KPoly<Rational> b) e) => new(e.a, e.b);
    public static implicit operator (KPoly<Rational> a, KPoly<Rational> b)(BGVCipher e) => (e.A, e.B);
}