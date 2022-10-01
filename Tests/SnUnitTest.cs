using FastGoat;
using Xunit;

namespace Tests;

public class SnUnitTest
{
    [Fact]
    public void Test1Tuple2Array()
    {
        Tuple2Array c1 = 1;
        Tuple2Array c2 = (2, 1);
        Tuple2Array c3 = (2, 3, 5);
        Tuple2Array c4 = (2, 3, 4, 5);
        Tuple2Array c5 = (8, 6, 2, 3, 5);
        Tuple2Array c6 = (2, 9, 3, 1, 5, 4);
        Tuple2Array c7 = (2, 8, 6, 11, 3, 5, 0);
        Tuple2Array c8 = (2, 8, 6, 11, 3, 10, 5, 8, 0);
        Tuple2Array c9 = (2, (1, 7), 5);

        Assert.Single(c1.Table);
        Assert.Equal(2, c2.Table.Length);
        Assert.Equal(3, c3.Table.Length);
        Assert.Equal(4, c4.Table.Length);
        Assert.Equal(5, c5.Table.Length);
        Assert.Equal(6, c6.Table.Length);
        Assert.Equal(7, c7.Table.Length);
        Assert.Empty(c8.Table);
        Assert.Empty(c9.Table);
    }
}