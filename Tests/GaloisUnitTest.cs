using System;
using System.Linq;
using FastGoat.UserGroup;
using Xunit;

namespace Tests;

public class GaloisUnitTest
{
    [Fact]
    public void Test1Ops()
    {
        foreach (var p in new[] { 2, 3, 5 })
        {
            for (int m = 1; m < 4; ++m)
            {
                var q = (int)Math.Pow(p, m); // q in 2,4,8,3,9,27,5,25,125
                var gf = FG.Galois(q);
                foreach (var e1 in gf)
                {
                    var ei = gf.Invert(e1);
                    Assert.Contains(ei, gf);
                    var e1ei = gf.Op(e1, ei);
                    Assert.Equal(gf.Neutral(), e1ei);
                    foreach (var e2 in gf)
                    {
                        Assert.Contains(gf.Op(e1, e2), gf);
                    }
                }
            }
        }
    }
}